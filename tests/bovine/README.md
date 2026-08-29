# Bovine imputation accuracy

Beagle's default `ne=100000` (effective population size) is calibrated for
human panels. Cattle breeds have an effective population size around 100,
and the imputation HMM's state-switch probability scales with
`0.04 * ne / nRefHaps` per cM — so at the default, the model switches
reference haplotypes far too readily and throws away the long shared
segments that make livestock imputation easy. Fixing `ne` is the single
largest accuracy lever we found; `preset=cattle` applies it.

```
rusty-beagle gt=chip.vcf.gz ref=panel.bref3 out=imputed preset=cattle nthreads=8
```

`preset=cattle` (alias `bovine`) currently expands to `ne=1000`. Explicitly
passed parameters always override the preset, and the run log records what
the preset filled in. Java Beagle users get the same benefit by passing
`ne=1000` directly — the preset is a convenience, not an algorithm change,
and without `preset=` rusty-beagle's output remains bit-identical to Java
Beagle 5.5.

## Evidence

### Real data: 500 dairy bulls (synbreedData, GPL-2)

`fetch_data.sh` downloads the public synbreedData cattle panel — 500 real
dairy bulls genotyped at 7,250 SNPs across all 29 autosomes, with pedigree —
from the CRAN GitHub mirror (md5-verified). `make_folds.py` builds 5
cross-validation folds that keep whole paternal half-sib families together
(random splits leak family haplotypes between reference and target and
overstate accuracy). Per fold, ~400 bulls form the reference (phased once
with Java Beagle, since reference-free phasing is not ported) and ~100
bulls are masked down to a low-density panel, imputed, and scored on the
withheld markers: squared correlation of the DS dosage with the true allele
count, computed per marker and averaged, plus best-guess concordance.

Mean masked-marker dosage r² (5 folds, sd ≈ 0.01):

| ne     | LD = every 5th marker | LD = every 10th marker |
|--------|----------------------|------------------------|
| 100000 (default) | 0.358 | 0.304 |
| 10000  | 0.707 | — |
| 5000   | 0.762 | — |
| 2000   | 0.772 | 0.496 |
| 1000   | 0.768 | 0.493 |
| 500    | 0.765 | 0.481 |
| 100    | 0.759 | — |

Concordance moves the same way (0.70 → 0.91 on the every-5th panel).
Everything else we tried — `err` from 0.0001 to 0.02, `cluster=0`,
whole-genome windows, `imp-nsteps=14`, `imp-segment=12`, re-phasing the
reference with tuned `ne` — shifted paired per-fold r² by at most 0.003.

### Simulated data: MacLeod et al. (2013) Holstein demography

`sim_cattle.py` simulates a cattle-like population with msprime using the
published Holstein Ne trajectory (Ne=90 today rising to 62,000 at 6,000
generations, as encoded by AlphaSimR's CATTLE species history; mu=2.5e-8,
rec=1e-8). 50 Mb chromosome, HD chip ascertained at 50k-chip density from
MAF>=0.05 variants, targets carrying every 10th HD marker, 250 target
animals, nested reference panels. Masked-marker dosage r²:

| ne     | 800 ref haps | 4,000 ref haps | 20,000 ref haps |
|--------|-------------|----------------|-----------------|
| 200    | 0.780 | 0.827 | 0.827 |
| 500    | **0.797** | 0.816 | 0.825 |
| 1000   | 0.789 | 0.828 | 0.831 |
| 2000   | 0.771 | 0.822 | 0.837 |
| 5000   | 0.681 | **0.829** | 0.834 |
| 20000  | 0.365 | 0.762 | **0.840** |
| 100000 (default) | 0.290 | 0.535 | 0.781 |

Two structural facts, consistent with the real data:

- the optimum is a broad plateau whose low side (ne 200–2000) is safe at
  every panel size, while the high side collapses for small panels;
- the penalty of the human default shrinks as the panel grows (because
  `ne/nRefHaps` falls), but even at 20,000 haplotypes it costs ~0.05 r².

ne=1000 is within ~0.01 of the measured peak in every condition, real and
simulated, which is why the preset uses it. At the 20,000-haplotype panel,
`err=0.01` and `imp-states=3200` each added <= 0.005 r² but showed no gain
on the real data, so the preset leaves them at Beagle defaults.

## Reproducing

```
# real data (needs R for the RData decode, Java Beagle for reference phasing)
tests/bovine/fetch_data.sh workdir/
python3 tests/bovine/make_folds.py --data-dir workdir --out-dir workdir/folds
for k in 0 1 2 3 4; do java -jar beagle.jar gt=workdir/folds/ref_$k.vcf.gz \
    out=workdir/folds/refph_$k nthreads=4 seed=1; done
rusty-beagle gt=workdir/folds/targ_0.vcf.gz ref=workdir/folds/refph_0.vcf.gz \
    out=imp0 preset=cattle nthreads=4
python3 tests/bovine/score.py --imputed imp0.vcf.gz \
    --truth workdir/folds/truth_0.npz --label "fold0 preset=cattle"

# simulated data
python3 tests/bovine/sim_cattle.py --out-dir simdir --length-mb 50 \
    --ref-sizes 400 2000 10000 --n-targ 250
rusty-beagle gt=simdir/simtarg.vcf.gz ref=simdir/simref_10000.vcf.gz \
    out=simimp preset=cattle nthreads=4
python3 tests/bovine/score.py --imputed simimp.vcf.gz --truth simdir/simtruth.npz \
    --label "sim preset=cattle"
```
