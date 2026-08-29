"""Scores an imputed VCF against withheld genotypes.

Reports, over the masked markers only:
  dosage r2   squared correlation of DS with the true allele count, computed
              per marker across target animals, averaged (the standard
              imputation-accuracy metric; insensitive to allele frequency
              inflation that concordance suffers from)
  concordance best-guess genotype (from GT) agreement with truth
binned by reference-panel allele frequency (MAF).
"""
import argparse, gzip, sys
import numpy as np

def load_imputed(path, want):
    """Returns {(chr,pos): (ds_row, gt_row)} for wanted markers, plus samples."""
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            p = line.rstrip("\n").split("\t")
            key = (int(p[0]), int(p[1]))
            if key not in want:
                continue
            fmt = p[8].split(":")
            i_ds = fmt.index("DS") if "DS" in fmt else None
            ds, gt = [], []
            for cell in p[9:]:
                sub = cell.split(":")
                a = sub[0].replace("|", "/").split("/")
                gt.append(int(a[0]) + int(a[1]))
                ds.append(float(sub[i_ds]) if i_ds is not None else float("nan"))
            out[key] = (np.array(ds), np.array(gt, dtype=np.int8))
    return out, samples

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imputed", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    t = np.load(a.truth, allow_pickle=True)
    truth, samples = t["truth"], list(t["samples"])
    keys = list(zip(t["chrom"].tolist(), t["pos"].tolist()))
    imputed, vcf_samples = load_imputed(a.imputed, set(keys))
    if vcf_samples != samples:
        order = [vcf_samples.index(s) for s in samples]
    else:
        order = None

    af = t["ref_af"]
    maf = np.minimum(af, 1 - af)
    bins = [(0.0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5001)]
    r2s = np.full(len(keys), np.nan)
    conc_num = np.zeros(len(keys)); conc_den = np.zeros(len(keys))
    for i, key in enumerate(keys):
        if key not in imputed:
            continue
        ds, gt = imputed[key]
        if order is not None:
            ds, gt = ds[order], gt[order]
        tr = truth[:, i].astype(np.float64)
        ok = tr >= 0
        if ok.sum() < 5:
            continue
        x, y = ds[ok], tr[ok]
        conc_num[i] = (gt[ok] == tr[ok]).sum(); conc_den[i] = ok.sum()
        if np.std(x) > 0 and np.std(y) > 0:
            r2s[i] = np.corrcoef(x, y)[0, 1] ** 2

    def fmt(sel):
        r = r2s[sel]; r = r[~np.isnan(r)]
        cn, cd = conc_num[sel].sum(), conc_den[sel].sum()
        return (f"r2={np.mean(r):.4f} (n={len(r)})  "
                f"conc={cn / max(cd, 1):.4f}") if len(r) else "n/a"

    sel_all = np.ones(len(keys), dtype=bool)
    print(f"[{a.label}] ALL     {fmt(sel_all)}")
    for lo, hi in bins:
        sel = (maf >= lo) & (maf < hi)
        print(f"[{a.label}] maf {lo:.2f}-{hi:.2f} ({sel.sum():4d}) {fmt(sel)}")
    # machine-readable summary line
    r = r2s[~np.isnan(r2s)]
    print(f"SUMMARY\t{a.label}\t{np.mean(r):.5f}\t{conc_num.sum()/max(conc_den.sum(),1):.5f}")

if __name__ == "__main__":
    main()
