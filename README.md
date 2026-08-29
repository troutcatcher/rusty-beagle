# rusty-beagle

A fast Rust port of [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html)
genotype phasing and imputation (upstream release `beagle.250227` /
`beagle.27Feb25.75f.jar`).

rusty-beagle produces **bit-identical output** to Java Beagle for
reference-based phasing and imputation — every phased GT, DS, GP/AP value and
every DR2/AF INFO field matches, byte for byte (only the `##filedate`/`##source`
header lines differ) — while running **1.4–2.9× faster** and using about
**2–3× less memory**.

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
panel (`ref=`, VCF format) and a phased *or unphased* target (`gt=`):

- gzip/BGZF VCF input, with parallel BGZF block decompression
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

Not (yet) ported: reference-free phasing — `ref=` is required, and `bref3`
reference files are not supported (use `unbref3` to convert to `.vcf.gz`).

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

All suites produce byte-identical output to `beagle.27Feb25.75f.jar`.
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
| 2,000 ref / 100 targ samples, 20k markers, 1 window | 7.7 s | 3.5 s | 2.2× |
| same, 12 windows (`window=4 overlap=2`) | 6.3 s | 2.2 s | 2.9× |
| 5,000 ref / 200 targ samples, 50k markers | 26.3 s | 15.5 s | 1.7× |
| peak RSS (5,000-sample run) | 1.32 GB | 0.64 GB | 2.1× less |

Phasing + imputation of an unphased target (same panels, ~5% of reference
markers genotyped in the target):

| workload | Java Beagle 5.5 | rusty-beagle | speedup |
|---|---|---|---|
| 2,000 ref / 100 targ samples, 20k markers | 12.8 s | 7.0 s | 1.8× |
| 5,000 ref / 200 targ samples, 50k markers | 43.4 s | 30.2 s | 1.4× |
| peak RSS (2,000-sample run) | 0.71 GB | 0.24 GB | 2.9× less |
| peak RSS (5,000-sample run) | 1.35 GB | 0.67 GB | 2.0× less |

Speed comes from the same parallel structure as Java (per-haplotype HMM,
per-sample phasing, per-cluster output, parallel input parsing) plus
Rust-side wins: no JVM warmup/GC, parallel BGZF inflation, flat
cache-friendly buffers, a bitset for allele matches, `Arc`-shared reference
haplotypes across phasing iterations, and zlib-rs for (de)compression.
Building with `RUSTFLAGS="-C target-cpu=native"` shaves off a further ~5–8%.
The phasing headroom is smaller than imputation's because most of the
phasing stage is spent in the `PhaseBaum2` forward/backward passes, whose
f32 evaluation order is pinned by bit-parity with Java.

## Build

```
cargo build --release          # target/release/rusty-beagle
cargo test                     # unit tests (java.util.Random port, formats)
```

## Validation harness

```
python3 tests/gen_test_data.py --out-dir /tmp/t1
bash tests/compare_beagle.sh /path/to/beagle.27Feb25.75f.jar /tmp/t1 nthreads=4
```

## License

GPL-3.0-or-later, the same license as Beagle. This project is a derivative
work of Beagle 5.5, Copyright (C) 2014–2024 Brian L. Browning
(`https://faculty.washington.edu/browning/beagle/beagle.html`). See
`docs/PORT_NOTES.md` for the source-to-source mapping.
