# rusty-beagle port notes

Port of Beagle 5.5 (source release `beagle.250227` / `beagle.27Feb25.75f.jar`)
genotype imputation to Rust. The reference Java source used for the port is the
Debian/Ubuntu packaging of the identical upstream release
(`beagle_250227.orig.tar.xz`, GPL-3+), which repackages
`https://faculty.washington.edu/browning/beagle/beagle.250227.zip`.

## Scope

v1 ports the **imputation** pipeline: `ref=` (phased reference VCF) +
`gt=` (phased target VCF) → imputed `out.vcf.gz`. This is exactly the code path
Java Beagle takes when every target genotype is phased and non-missing
(`Main.phaseAndImpute` skips phasing whenever `fpd.targGT().isPhased()`), so
the two programs are directly comparable. Unphased/missing target genotypes
require the phasing stage (Java `phase` package), which is not ported yet;
rusty-beagle exits with a clear error in that case.

The goal is *bit-identical output* to Java Beagle (up to the `##filedate`
header line), achieved by porting the algorithms operation-for-operation,
including:

- `java.util.Random` (48-bit LCG) — `src/javautil.rs`
- `java.util.PriorityQueue` binary-heap sift-up/down order (ties in
  `CompHapSegment` ordering are resolved by heap layout, which affects
  composite-haplotype construction) — `src/javautil.rs`
- `blbutil.Utilities.shuffle` partial Fisher–Yates
- float (f32) arithmetic in the HMM forward/backward pass in the same
  evaluation order as Java
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
| `vcf.RefTargSlidingWindow` | `windows.rs` | window/overlap/splice logic, `MarkerIndices` |
| `vcf.GeneticMap`/`PlinkGenMap`/`PositionMap` | `genmap.rs` | PLINK map interpolation w/ 5cM end rule; default pos*1e-6 |
| `imp.ImpData` + `imp.HaplotypeCoder` | `impdata.rs` | marker clusters, per-cluster allele-sequence coding (both `codeSeq` and the seq-coded composition path) |
| `imp.CodedSteps` + `imp.ImpIbs` | `impibs.rs` | PBWT-like partition refinement per step; `Random(seed + parent[0])` subsets |
| `imp.ImpStates` + `beagleutil.CompHapSegment` | `impstates.rs` | composite reference haplotypes via Java-ordered priority queue |
| `imp.ImpLSBaum` + `imp.StateProbsFactory` | `hmm.rs` | Li–Stephens fwd/bwd in f32; state-prob sparsification threshold `min(0.005, 0.9999/nStates)` |
| `imp.ImputedVcfWriter` + `imp.RefHapHash` + `imp.ImputedRecBuilder` | `impout.rs` | per-cluster output, allele-seq hashing (`Random(start)`), DS/AP/GP/DR2/AF formatting |
| `main.WindowWriter` + `blbutil.BGZIPOutputStream` | `bgzf.rs`, `impout.rs` | BGZF blocks compressed in parallel, final empty BGZF block |

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

## Validation

`tests/` contains a harness that generates synthetic phased ref/target panels,
masks markers from the target, runs both Java Beagle (`beagle.jar`, from the
Debian package of the same release) and rusty-beagle, and diffs the gunzipped
output VCFs (ignoring `##filedate`). See README for results.
