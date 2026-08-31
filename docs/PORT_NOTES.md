# rusty-beagle port notes

Port of Beagle 5.5 (source release `beagle.250227` / `beagle.27Feb25.75f.jar`)
genotype imputation to Rust. The reference Java source used for the port is the
Debian/Ubuntu packaging of the identical upstream release
(`beagle_250227.orig.tar.xz`, GPL-3+), which repackages
`https://faculty.washington.edu/browning/beagle/beagle.250227.zip`.

## Scope

rusty-beagle ports the **phasing + imputation** pipeline: `ref=` (phased
reference VCF or bref3) + `gt=` (phased or unphased target VCF) →
phased/imputed `out.vcf.gz`. A fully phased, non-missing target takes the
fast path (`Main.phaseAndImpute` skips phasing whenever
`fpd.targGT().isPhased()`); otherwise the target is phased against the
reference panel first, exactly as in the Java `phase` package. `impute=false`
restricts windows to target markers and prints phased target genotypes only.
Reference-free phasing (no `ref=`) is not ported; rusty-beagle exits with a
clear error in that case, as does old-format `.bref` (non-v3) input, matching
Java Beagle's own rejection of that format.

The goal is *bit-identical output* to Java Beagle (up to the `##filedate`
header line), achieved by porting the algorithms operation-for-operation,
including:

- `java.util.Random` (48-bit LCG) — `src/javautil.rs`
- `java.util.PriorityQueue` binary-heap sift-up/down order (ties in
  `CompHapSegment` ordering are resolved by heap layout, which affects
  composite-haplotype construction) — `src/javautil.rs`
- `blbutil.Utilities.shuffle` partial Fisher–Yates
- float (f32) arithmetic in the HMM forward/backward passes in the same
  evaluation order as Java, including Java's `float += double` compound
  assignment (add in f64, then narrow)
- `blbutil.BitArray`'s XOR range hash, which depends on the absolute bit
  alignment of alleles within the backing u64 words — `src/bits.rs`
- `DecimalFormat("#.##")`-style tables for DS/GP values, DR2 (2dp) and
  AF (4dp) rounding
- `Math.rint` = round-half-to-even for dose indexing

## Pipeline map (Java → Rust)

| Java | Rust | Notes |
|---|---|---|
| `main.Main`, `main.Par` | `main.rs`, `par.rs` | key=value CLI, same defaults |
| `blbutil.InputIt`/`BGZipIt` | `vcfio.rs` | gzip/BGZF text input (flate2 MultiGzDecoder) |
| `vcf.VcfHeader`, `vcf.Samples` | `vcfio.rs` | diploid/haploid from first data line |
| `vcf.Marker`+`MarkerParser` | `marker.rs` | stores ID + REF/ALT + INFO/END; equality = (chrom,pos,alleles,END) |
| `vcf.VcfIt` (`TO_LOWMEM_GT_REC`) | `vcfio.rs` | target GT records; phased flag per record |
| `vcf.RefIt` + `bref.SeqCoder3` | `refpanel.rs` | allele-coded records + sequence-coding simulation; the seq-coded *block boundaries* (flush on chrom change / EOF / >maxNSeq sequences) are semantically relevant: they force marker-cluster splits in `ImpData.targBlockEnd` |
| `bref.Bref3Header` + `bref.Bref3Reader` + `bref.Bref3It` | `bref3.rs` | binary reference (`.bref3`) input; decodes the file's own stored block/sequence-coding structure directly (no VCF round-trip, no re-run of `SeqCoder3`) — see notes below |
| `vcf.RefTargSlidingWindow` | `windows.rs` | window/overlap/splice logic, `MarkerIndices`; dispatches `ref=` to `refpanel::RefReader` or `bref3::Bref3RefReader` behind a shared `RefSource` trait, by file extension, exactly like Java's `.bref3` vs. VCF dispatch |
| `vcf.GeneticMap`/`PlinkGenMap`/`PositionMap` | `genmap.rs` | PLINK map interpolation w/ 5cM end rule; default pos*1e-6 |
| `imp.ImpData` + `imp.HaplotypeCoder` | `impdata.rs` | marker clusters, per-cluster allele-sequence coding (both `codeSeq` and the seq-coded composition path) |
| `imp.CodedSteps` + `imp.ImpIbs` | `impibs.rs` | PBWT-like partition refinement per step; `Random(seed + parent[0])` subsets |
| `imp.ImpStates` + `beagleutil.CompHapSegment` | `impstates.rs` | composite reference haplotypes via Java-ordered priority queue |
| `imp.ImpLSBaum` + `imp.StateProbsFactory` | `hmm.rs` | Li–Stephens fwd/bwd in f32; state-prob sparsification threshold `min(0.005, 0.9999/nStates)` |
| `imp.ImputedVcfWriter` + `imp.RefHapHash` + `imp.ImputedRecBuilder` | `impout.rs` | per-cluster output, allele-seq hashing (`Random(start)`), DS/AP/GP/DR2/AF formatting |
| `main.WindowWriter` + `blbutil.BGZIPOutputStream` | `bgzf.rs`, `impout.rs` | BGZF blocks compressed in parallel, final empty BGZF block |
| `blbutil.BitArray` | `bits.rs` | u64-word bit list; range hash/equality/copy/swap with Java's exact word masking |
| `phase.FixedPhaseData` + `phase.PhaseData` | `phasedata.rs` | spliced target alleles, rare-allele carriers, hi-freq marker set, stage-1 restriction, per-iteration seeds, `pMismatch`/`recombIntensity` updates |
| `phase.SamplePhase` | `phasedata.rs` | per-sample genotype clusters (missing / masked-het / hom / phased-het / unphased-het), trailing-unphased-het masking |
| `vcf.XRefGT` | `xref.rs` | hap-major bit-packed haplotypes; the combined target+reference view shares reference haps via `Arc` instead of copying them each iteration |
| `phase.CodedSteps` | `codedsteps.rs` | per-step unique allele-sequence indexing (first-seen order) |
| `phase.Ibs2`, `phase.Ibs2Sets`, `phase.Ibs2Markers` | `ibs2.rs` | IBS2 segment detection used to exclude relatives from phasing states |
| `phase.PbwtPhaser`/`PbwtRecPhaser`/`RevPbwtPhaser`/`FwdPbwtPhaser` | `initphase.rs`, `pbwt.rs` | PBWT initial phasing (forward + reverse passes with `Random` tie-breaks) |
| `phase.PbwtPhaseIbs` + `phase.LowFreqPbwtPhaseIbs` | `phaseibs.rs` | batched PBWT IBS candidate selection for stage 1/2, with buffer regions and backoff |
| `phase.BasicPhaseStates` + `phase.PhaseBaum2` | `phasebaum.rs` | composite-haplotype state lists and the three-track phasing HMM (het phasing, missing/masked imputation, hap swap tracking) |
| `phase.HmmParamData` + `phase.ParamEstimates` | `phasebaum.rs` | EM estimation of `pMismatch`/`recombIntensity` (estimates are sorted before summing, which makes results thread-count-stable) |
| `phase.LowFreqPhaseStates` + `phase.HmmStateProbs` + `phase.Stage2Baum`/`Stage2Haps` | `stage2.rs` | second-stage phasing of rare variants at stage-1 marker anchors |

## Behavioural notes discovered while porting

- `ImpData.MIN_CM_DIST = 1e-7` is applied to *cumulative* target/ref genetic
  positions.
- `pRecomb = -expm1(-0.04*ne/nRefHaps * dPos)` in f32.
- `err` default is the Li–Stephens/`Marchini et al.` approximation
  `θ/(2(θ+H))`, `θ = 1/(ln(H)+0.5)`, H = nRefHaps+nTargHaps; per-cluster error
  = `err * clusterSize`, capped at 0.5.
- Window splice points: `prevSplice = overlapEnd>>1`,
  `nextSplice = (nMarkers+overlapStart)>>>1`, and their target-index versions
  via binary search.
- The `window` reader drops target markers that are absent from the reference
  panel (same position+alleles+END required).
- `SeqCoder3.defaultMaxNSeq(n) = floor(2^(2*log10(n)+1))` clamped to 65535;
  a record is seq-coded iff `nAlleles <= min(maxNSeq,255)` and major-allele
  count `<= floor(0.995*nHaps - 1)`.
- `ImpIbs`: haplotypes per step = `imp_states / round(imp_segment/imp_step)`;
  steps merged per segment = `imp_nsteps`.
- `ImpStates`: when an IBS hap is added and the queue is full, the head is
  recycled; the segment switch point is `stepStart((headLastStep+curStep)/2)`.
- `StateProbsFactory` stores P(state) at marker m and m+1 for interpolation.
- Output DS/GP values print from a 201-entry `#.##` table; DR2 prints with
  exactly 2 decimals, AF with 4 (`0.0000` format, HALF_EVEN); hom-ref samples
  short-circuit to a shared `0|0:0[,0...]` string. DR2 formula per ALT allele
  uses f32 accumulation over haplotype doses in sample order.
- The output INFO field is `DR2=..;AF=..[;END=..][;IMP]`, FORMAT `GT:DS[:AP1:AP2][:GP]`.
- Monomorphic output records (`nAlleles==1`) print INFO `IMP`/empty and GT only.

Phasing-specific notes:

- Phased output depends on (`seed`, `nthreads`), in Java and in the port:
  `nthreads` enters through `PbwtPhaser`'s hi-freq window partitioning
  (`advanceCM = max(2.0, totalCM/nThreads)`), `PbwtIbsData.stepsPerBatch`,
  and `CodedSteps`' batch splitting. For a fixed pair the run is fully
  deterministic.
- The EM early-stop `hpd.sumSwitchProbs() < maxSum` in
  `PhaseLS.getParamEstimates` is dead code: `addEstimationData` drains the
  accumulator before the next loop check, so every worker processes samples
  until the shared counter is exhausted. `ParamEstimates` sorts its
  accumulated lists before summing, which is why thread scheduling does not
  affect the estimates. Both quirks are preserved.
- Burn-in ends early when the per-iteration hap-swap rate drops to ≤ 0.01
  (`PhaseData.advanceToFirstPhasingIt`).
- Genotype clusters within a sample split at 0.005 cM or 255 genotypes;
  trailing unphased hets within 3000 bp of the window end are masked and
  re-imputed (`SamplePhase`).
- Stage 2 (rare variants) runs only when some target markers are excluded
  from stage 1; its per-sample RNG is `rand.setSeed(itSeed + sample)`.
- Windows with no ungenotyped reference markers (and `impute=false` runs)
  print phased target genotypes without dosage fields, via
  `WindowWriter.printPhased`.

bref3-specific notes:

- A bref3 file is a sequence of self-describing blocks (record count, chrom
  name, a shared hap→seq map, then per-marker data) ending in a zero-record
  sentinel, followed by a random-access index that neither Java's `Bref3It`
  nor rusty-beagle ever reads — both decode strictly sequentially from the
  start of the file.
- Each on-disk block becomes one `RefAlleles::SeqCoded` block id in
  `RefRec`, reproducing Java's `rec.map(0) != hapToSeq` object-identity
  check: every marker decoded from one block-header read shares one id
  (one shared `hap2seq` `Arc`), exactly as they share one `hapToSeq`
  reference in Java. rusty-beagle does not re-run `SeqCoder3` on bref3
  input — the file's stored block boundaries are used verbatim.
- A marker is stored either sequence-coded (a shared per-block hap→seq map
  plus a per-record seq→allele byte array) or, when it didn't compress well
  at conversion time, allele-coded (a length-prefixed haplotype list per
  non-major allele, with a `-1` length sentinel marking the major allele's
  slot).
- Because those block boundaries are fixed once, at `bref3` conversion
  time, by whatever `SeqCoder3` state the *conversion* run happened to be
  in, they are not guaranteed to match the boundaries a fresh read of the
  same data from VCF would independently choose. Since block boundaries
  force marker-cluster splits (`ImpData.targBlockEnd`), imputing from a
  pre-converted `.bref3` file can therefore give different — but equally
  valid — output than imputing from the VCF it was converted from. This is
  a property of the Java algorithm, reproduced faithfully rather than
  worked around; verified empirically (Java itself gives different output
  for the two inputs on a multi-chromosome test panel), so the bref3
  parity suites compare rusty-beagle against Java reading the *same*
  `.bref3` file rather than against a VCF baseline.
- `excludesamples=`/`excludemarkers=`/`chrom=` are applied when the bref3
  file is *read*, independent of how it was converted: sample exclusion
  remaps haplotype indices (`Bref3Header.filteredHapIndices` /
  `invfilteredHapIndices`) while decoding each block's hap→seq map and each
  allele-coded record's per-allele haplotype lists; excluded samples' data
  is skipped, not decoded into the panel.
- Sample IDs, marker IDs, chromosome names, and REF/ALT alleles are stored
  as Java's "modified UTF-8" (`DataInput.readUTF`): standard UTF-8 except
  NUL is coded as the two bytes `0xC0,0x80`, and characters outside the
  BMP are coded as a pair of independently 3-byte-encoded UTF-16
  surrogates (CESU-8) rather than one 4-byte sequence. `bref3.rs` decodes
  this explicitly rather than assuming plain UTF-8, covered by unit tests
  since no real-world genomic data exercises those code paths.
- SNV alleles (a REF/ALT set drawn from A/C/G/T) are coded compactly as one
  byte: a permutation index into the 24 orderings of "A","C","G","T" plus
  an allele count, rather than as explicit strings. The 24-entry
  permutation table is cross-checked by unit test against a direct port of
  `Bref3Reader.snvPerms()`'s recursive generator.

## Validation

`tests/` contains a harness that generates synthetic phased ref/target panels,
masks markers from the target, runs both Java Beagle (`beagle.jar`, from the
Debian package of the same release) and rusty-beagle, and diffs the gunzipped
output VCFs (ignoring `##filedate`). `tests/compare_beagle_bref3.sh` is the
same idea for `.bref3` reference input: it converts `ref.vcf.gz` with the
separate `bref3` conversion tool, then runs both programs with
`ref=<the bref3 file>` (not against the VCF baseline — see the bref3-specific
notes above for why that comparison isn't meaningful). See README for
results.

## License and provenance

Because this document describes a class-by-class translation, it is also the
record of *why* rusty-beagle must be GPL-licensed: the port is a derivative
work of Beagle 5.5, Copyright (C) 2014-2024 Brian L. Browning, which is
distributed under the GNU General Public License, version 3 or later. Section
5 of that license requires modified versions to be released under the same
terms, with prominent notices that they are modified, so:

- `LICENSE` holds the complete GPL v3 text (GPL section 4 requires shipping it);
- `NOTICE` holds the copyright, the upstream attribution, the statement that
  the whole of this work is modified relative to Beagle 5.5, and the
  third-party crate licenses;
- every file under `src/` and `tests/` carries an SPDX tag and the GPL notice,
  and names the upstream class it corresponds to;
- `rusty-beagle` prints the copyright and warranty notice at startup, writes
  them into each run's `.log`, and prints the full notice for
  `rusty-beagle license`.

The mapping tables above satisfy the "prominent notices stating that you
modified it" requirement in substance as well as in form: each entry states
which upstream class a Rust module was derived from, and the behavioural notes
state where the port's behaviour was deliberately changed.
