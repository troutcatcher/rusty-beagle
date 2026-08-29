# rusty-beagle

A fast Rust port of [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html)
genotype imputation (upstream release `beagle.250227` / `beagle.27Feb25.75f.jar`).

rusty-beagle produces **bit-identical output** to Java Beagle for
reference-based imputation of a phased target — every GT, DS, GP/AP value and
every DR2/AF INFO field matches, byte for byte (only the `##filedate`/`##source`
header lines differ) — while running about **2–3× faster** and using about
**2.5× less memory**.

```
rusty-beagle gt=target.vcf.gz ref=reference.vcf.gz out=imputed map=plink.chr20.map nthreads=8
```

The command line is Beagle's: the same `key=value` parameters with the same
defaults (`window`, `overlap`, `ne`, `err`, `imp-states`, `imp-segment`,
`imp-step`, `imp-nsteps`, `cluster`, `ap`, `gp`, `chrom`, `excludesamples`,
`excludemarkers`, `map`, `seed`, `nthreads`, ...), so it drops into existing
pipelines.

## Scope

Ported: the full **imputation** pipeline for a phased reference panel
(`ref=`, VCF format) and a phased target (`gt=`):

- gzip/BGZF VCF input, with parallel BGZF block decompression
- marker windows, window overlap and splice logic
- PLINK genetic-map interpolation
- Beagle's low-memory allele-coded / sequence-coded reference representation
  (`SeqCoder3`), whose block boundaries also determine marker clustering
- IBS-based composite reference haplotypes (`ImpIbs`/`ImpStates`)
- the Li–Stephens forward–backward HMM (`ImpLSBaum`)
- linear interpolation of state probabilities at ungenotyped markers
- DS/AP/GP dosages, DR2/AF/IMP INFO fields, BGZF-compressed output

Not (yet) ported: **phasing**. Java Beagle phases an unphased target before
imputing; rusty-beagle requires the target to be phased and non-missing, and
exits with a clear message otherwise (phase once with Java Beagle
`impute=false`, then impute with rusty-beagle). When every target genotype is
phased, Java Beagle skips its phasing stage entirely, which is what makes the
two programs directly comparable. `bref3` reference files are also not
supported yet — use `unbref3` to convert to `.vcf.gz`.

## Correctness

`tests/compare_beagle.sh` runs Java Beagle and rusty-beagle on the same input
and diffs the decompressed output VCFs. The committed harness generates
synthetic panels (`tests/gen_test_data.py`) covering:

| suite | features exercised |
|---|---|
| t1 | defaults, 190 ref / 10 targ samples, 1k markers |
| t2 | multi-window (`window=4 overlap=2`), multiallelic markers, haploid samples, INFO/END, `gp=true` |
| t3 | 1500 ref samples, `imp-states=64` (composite-haplotype recycling), `ap=true` |
| t4 | PLINK map, `chrom=22:300000-5500000`, `excludesamples`/`excludemarkers`, `ne`/`err`/`cluster` overrides |
| t5 | two chromosomes in one VCF |

All suites produce byte-identical output to `beagle.27Feb25.75f.jar`.
Bit-parity is achieved by porting the algorithms operation-for-operation,
including `java.util.Random`, `java.util.PriorityQueue`'s heap layout (its
tie-breaking affects composite-haplotype construction), f32 arithmetic order
in the HMM, and `DecimalFormat`/`Math.rint` rounding (see
`docs/PORT_NOTES.md`).

## Performance

Measured on 4 cores (`nthreads=4`), synthetic panels, wall time including
JVM/process startup; outputs verified identical in every run:

| workload | Java Beagle 5.5 | rusty-beagle | speedup |
|---|---|---|---|
| 2,000 ref / 100 targ samples, 20k markers, 1 window | 7.7 s | 3.5 s | 2.2× |
| same, 12 windows (`window=4 overlap=2`) | 6.3 s | 2.2 s | 2.9× |
| 5,000 ref / 200 targ samples, 50k markers | 26.3 s | 15.5 s | 1.7× |
| peak RSS (5,000-sample run) | 1.32 GB | 0.64 GB | 2.1× less |

Speed comes from the same parallel structure as Java (per-haplotype HMM,
per-cluster output, parallel input parsing) plus Rust-side wins: no JVM
warmup/GC, parallel BGZF inflation, flat cache-friendly buffers, a bitset for
allele matches, and zlib-rs for (de)compression. Building with
`RUSTFLAGS="-C target-cpu=native"` shaves off a further ~5–8%.

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
