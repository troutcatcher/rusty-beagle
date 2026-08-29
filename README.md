# rusty-beagle

A fast Rust port of [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html)
genotype phasing and imputation (upstream release `beagle.250227` /
`beagle.27Feb25.75f.jar`).

rusty-beagle produces **bit-identical output** to Java Beagle for
reference-based phasing and imputation — every phased GT, DS, GP/AP value and
every DR2/AF INFO field matches, byte for byte (only the `##filedate`/`##source`
header lines differ) — while running **2–6× faster** and using **2–8× less
memory**, depending on the shape of the panel.

```
rusty-beagle gt=target.vcf.gz ref=reference.vcf.gz out=imputed map=plink.chr20.map nthreads=8
```

The command line is Beagle's: the same `key=value` parameters with the same
defaults (`window`, `overlap`, `ne`, `err`, `imp-states`, `imp-segment`,
`imp-step`, `imp-nsteps`, `cluster`, `ap`, `gp`, `chrom`, `excludesamples`,
`excludemarkers`, `map`, `seed`, `nthreads`, ...), so it drops into existing
pipelines.

## Scope

Ported: the full **phasing + imputation** pipeline for a phased reference
panel (`ref=`, VCF or bref3 format) and a phased *or unphased* target
(`gt=`):

- gzip/BGZF VCF input, with parallel BGZF block decompression
- binary reference (`.bref3`) input, decoded directly with no VCF
  round-trip — see below
- marker windows, window overlap and splice logic
- PLINK genetic-map interpolation
- Beagle's low-memory allele-coded / sequence-coded reference representation
  (`SeqCoder3`), whose block boundaries also determine marker clustering
- the phasing stage (`phase` package): PBWT-based initial phasing, burn-in +
  phasing iterations of the three-track Li–Stephens HMM (`PhaseBaum2`) with
  EM estimation of the mismatch/recombination parameters, IBS2 constraint
  segments, and second-stage phasing of rare variants; `impute=false`
  produces phased-target-only output
- IBS-based composite reference haplotypes (`ImpIbs`/`ImpStates`)
- the imputation Li–Stephens forward–backward HMM (`ImpLSBaum`)
- linear interpolation of state probabilities at ungenotyped markers
- DS/AP/GP dosages, DR2/AF/IMP INFO fields, BGZF-compressed output

Like Java Beagle, phased output depends on `seed` *and* `nthreads` (Beagle
partitions PBWT windows/batches by thread count); rusty-beagle reproduces
Java's output exactly for any given (`seed`, `nthreads`) pair.

`ref=some.bref3` is a drop-in alternative to `ref=some.vcf.gz`: rusty-beagle
decodes the file's own block/sequence-coding structure directly instead of
re-deriving it from VCF text, which is both faster and lower-memory than the
VCF path (bref3 is already compact, so there's no decompression or text
parsing to do). Note that a bref3 file's block boundaries are fixed at
*conversion* time and are semantically relevant to marker clustering (see
`docs/PORT_NOTES.md`), so imputing from a pre-converted bref3 file can give
different (equally valid) output than imputing from the original VCF the
bref3 was built from — this is inherent Beagle behavior, not a rusty-beagle
quirk, and rusty-beagle matches Java Beagle exactly for either input.

Not (yet) ported: reference-free phasing — `ref=` is required. Old-format
`.bref` (not `.bref3`) files are not supported, matching Java Beagle, which
also rejects them.

## Species presets: bovine accuracy

Beagle's `ne=100000` default is calibrated for human panels; cattle have an
effective population size near 100, and the imputation HMM's switch rate
scales with `0.04*ne/nRefHaps`. `preset=cattle` (a rusty-beagle extension;
explicit parameters always override it, and without it output stays
bit-identical to Java Beagle) applies `ne=1000`, which sits on the accuracy
plateau at every panel size we measured:

- 500 real dairy bulls (public synbreedData panel, sire-family-aware 5-fold
  CV, LD-chip targets imputed to the 7,250-marker panel): masked-marker
  dosage r² **0.36 → 0.77**, concordance 0.70 → 0.91.
- Coalescent panels under the published MacLeod et al. (2013) Holstein
  demography: r² +0.50 at an 800-haplotype reference, +0.05 at 20,000
  haplotypes.

Java Beagle users get the same gain by passing `ne=1000` explicitly.

On top of that, `ensemble=K` (also a rusty-beagle extension, default 1 =
bit-identical) runs the seeded-stochastic phasing/imputation K times
internally, aligns the replicates' haplotypes, and averages the allele
probabilities. At K=5 it adds a further +0.02 to +0.06 dosage r² on every
bovine dataset tested (largest on sparse chips and small references, where
phasing noise dominates), at ~K-fold run time. Unlike the preset this is a
genuine algorithm extension: single-run Java Beagle cannot reproduce it,
though averaging K separate Java runs' DS fields gets most of the way.

Full methodology, sweep tables and reproduction commands:
`tests/bovine/README.md`.

## Correctness

`tests/compare_beagle.sh` runs Java Beagle and rusty-beagle on the same input
and diffs the decompressed output VCFs. The committed harness generates
synthetic panels (`tests/gen_test_data.py`) covering:

| suite | features exercised |
|---|---|
| t1 | defaults, 190 ref / 10 targ samples, 1k markers (phased target) |
| t2 | multi-window (`window=4 overlap=2`), multiallelic markers, haploid samples, INFO/END, `gp=true` |
| t3 | 1500 ref samples, `imp-states=64` (composite-haplotype recycling), `ap=true` |
| t4 | PLINK map, `chrom=22:300000-5500000`, `excludesamples`/`excludemarkers`, `ne`/`err`/`cluster` overrides |
| t5 | two chromosomes in one VCF |

Phasing suites (unphased and partially phased targets):

| suite | features exercised |
|---|---|
| ph1 | unphased target with haploid samples and multiallelic markers; `nthreads=1` and `nthreads=4` each match the corresponding Java run |
| ph2 | determinism: three repeated runs with the same `seed`, all identical |
| ph3 | multi-window phasing (`window=4 overlap=2`), missing genotypes, multiallelics, `gp=true` |
| ph4 | 600 target samples vs 200 reference samples; second-stage (rare-variant) phasing |
| ph5 | phasing with PLINK map, `chrom`, excludes, `ne`/`err`/`cluster` overrides, `seed` |
| ph6 | haploid samples mixed with unphased diploid samples |
| ph7 | two chromosomes, multi-window phasing |
| — | `impute=false` (phase-only output), single- and multi-window |
| — | `.bref3` reference input (both default and forced-small-block encodings), including with `excludesamples`/`excludemarkers`/`chrom` and with phasing |

All suites produce byte-identical output to `beagle.27Feb25.75f.jar`. Since a
bref3 file's own block boundaries can make its output legitimately differ
from imputing the same data straight from VCF (see Scope above), the bref3
suites compare rusty-beagle against Java Beagle reading the *same* `.bref3`
file, generated with the `bref3` conversion tool from the same release.
Bit-parity is achieved by porting the algorithms operation-for-operation,
including `java.util.Random`, `java.util.PriorityQueue`'s heap layout (its
tie-breaking affects composite-haplotype construction), f32 arithmetic order
in the HMMs, Java `BitSet`-layout-dependent haplotype hashes, and
`DecimalFormat`/`Math.rint` rounding (see `docs/PORT_NOTES.md`).

## Performance

Measured on 4 cores (`nthreads=4`), synthetic panels, wall time including
JVM/process startup; outputs verified identical in every run.

Imputation of an already-phased target:

| workload | Java Beagle 5.5 | rusty-beagle | speedup |
|---|---|---|---|
| 2,000 ref / 100 targ samples, 20k markers, 1 window | 8.5 s | 2.4 s | 3.5× |
| same, 12 windows (`window=4 overlap=2`) | 5.7 s | 1.7 s | 3.4× |
| 5,000 ref / 200 targ samples, 50k markers | 27.4 s | 11.0 s | 2.5× |
| peak RSS (5,000-sample run) | 1.70 GB | 0.66 GB | 2.6× less |

Phasing + imputation of an unphased target (same panels, ~5% of reference
markers genotyped in the target):

| workload | Java Beagle 5.5 | rusty-beagle | speedup |
|---|---|---|---|
| 2,000 ref / 100 targ samples, 20k markers | 10.2 s | 3.6 s | 2.8× |
| 5,000 ref / 200 targ samples, 50k markers | 30.2 s | 14.2 s | 2.1× |
| peak RSS (2,000-sample run) | 0.93 GB | 0.22 GB | 4.2× less |
| peak RSS (5,000-sample run) | 1.60 GB | 0.64 GB | 2.5× less |

Imputation of a phased target from a `.bref3` reference panel (same panels,
converted with the `bref3` tool at its default block size) is faster still,
since there is no gzip decompression or VCF text parsing to do:

| workload | Java Beagle 5.5 | rusty-beagle | speedup |
|---|---|---|---|
| 2,000 ref / 100 targ samples, 20k markers | 7.0 s | 1.7 s | 4.1× |
| 5,000 ref / 200 targ samples, 50k markers | 22.2 s | 7.8 s | 2.8× |
| peak RSS (5,000-sample run) | 1.35 GB | 0.62 GB | 2.2× less |

A livestock-shaped panel — many animals, one chromosome, few SNPs — where
the reference is a `.bref3` file of 200,000 animals genotyped at 1,300 SNPs
across a 158 Mb chromosome, and the target cohort carries every 10th SNP:

| target cohort | Java Beagle 5.5 | rusty-beagle | speedup |
|---|---|---|---|
| 500 animals, imputation | 9.0 s | 1.6 s | 5.6× |
| 500 animals, phasing + imputation | 14.3 s | 4.2 s | 3.4× |
| 5,000 animals, imputation | 32.5 s | 8.1 s | 4.0× |
| 5,000 animals, phasing + imputation | 40.1 s | 10.1 s | 4.0× |
| 20,000 animals, imputation | 125.3 s | 32.6 s | 3.8× |
| 20,000 animals, phasing + imputation | 121.1 s | 30.6 s | 4.0× |
| peak RSS (500-animal run) | 2.18 GB | 0.26 GB | 8.3× less |
| peak RSS (20,000-animal run) | 12.60 GB | 4.59 GB | 2.7× less |

Peak memory in that regime is dominated by the retained per-haplotype HMM
state probabilities, which both programs hold for every target haplotype at
once, so it grows with the target cohort rather than with the reference
panel.

Speed comes from the same parallel structure as Java (per-haplotype HMM,
per-sample phasing, per-cluster output, parallel input parsing) plus
Rust-side wins: no JVM warmup/GC, parallel BGZF inflation, and zlib-rs for
(de)compression. Most of the rest is layout and allocation work in the hot
loops rather than anything algorithmic, since the arithmetic order is pinned
by bit-parity with Java:

- VCF lines are read in batches into one reused buffer, and line reading runs
  on its own thread so gunzip overlaps with parsing and sequence-coding
- composite-haplotype segment switches are precomputed into a cluster-sorted
  transition list, instead of rescanning all states at every cluster
- the per-cluster allele-match bitset is built one 64-bit word at a time with
  a branchless shift, and reads a precomputed match bitset (~1 KB, L1-sized)
  rather than gathering from the full per-haplotype code array
- per-haplotype ALT-allele lists and state probabilities use flat CSR buffers
  instead of `Vec<Vec<_>>` rebuilt per cluster, with the retained states
  interleaved so the values used together are read together
- the backward pass records which states survive sparsification while their
  probabilities are still in hand, so the retained rows are allocated once at
  exactly their final size
- a state's reference haplotype resolves to its table index through a direct
  array rather than a binary search, and the table's haplotypes are deduped
  with a bitset instead of sorting the full multiset (millions of entries per
  cluster, tens of thousands distinct)

Phasing repeats its reference-panel work on every iteration -- re-coding
allele sequences per step and rebuilding the PBWT over every haplotype -- so
against a large reference that work, not the per-sample HMM, is what it
costs. Both are addressed above (the integer hasher and the narrowed step
tables); what remains is the `PhaseBaum2` forward/backward passes, which are
float arithmetic in an evaluation order fixed by bit-parity.

`RUSTFLAGS="-C target-cpu=native"` is not recommended: on the benchmark
machine it came out marginally *slower* than the portable build once the hot
loops were branchless (output is identical either way). Measure before
adopting it.

## Build

```
cargo build --release          # target/release/rusty-beagle
cargo test                     # unit tests (java.util.Random port, formats)
```

## Validation harness

```
python3 tests/gen_test_data.py --out-dir /tmp/t1
bash tests/compare_beagle.sh /path/to/beagle.27Feb25.75f.jar /tmp/t1 nthreads=4

# large-cohort panel (200k animals, one chromosome), then convert to bref3
python3 tests/gen_large_cohort.py 200000 5000 1300 158000000 /tmp/big
java -jar /path/to/bref3.27Feb25.75f.jar /tmp/big/ref.vcf > /tmp/big/ref.bref3

# .bref3 reference input (needs the separate bref3 conversion tool)
bash tests/compare_beagle_bref3.sh /path/to/beagle.27Feb25.75f.jar /path/to/bref3.27Feb25.75f.jar /tmp/t1 nthreads=4
```

## License

GPL-3.0-or-later, the same license as Beagle. This project is a derivative
work of Beagle 5.5, Copyright (C) 2014–2024 Brian L. Browning
(`https://faculty.washington.edu/browning/beagle/beagle.html`). See
`docs/PORT_NOTES.md` for the source-to-source mapping.
