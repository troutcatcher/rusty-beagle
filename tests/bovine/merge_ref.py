"""Appends imputed-cohort haplotypes to a phased reference VCF (two-pass
cohort imputation): both inputs must cover the same markers in order."""
import gzip, sys
ref_path, imp_path, out_path = sys.argv[1:4]
def rows(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                yield line.rstrip("\n")
r, i = rows(ref_path), rows(imp_path)
with gzip.open(out_path, "wt") as out:
    out.write("##fileformat=VCFv4.2\n")
    for a, b in zip(r, i):
        pa, pb = a.split("\t"), b.split("\t")
        if pa[0] == "#CHROM":
            out.write("\t".join(pa + pb[9:]) + "\n")
            continue
        assert (pa[0], pa[1]) == (pb[0], pb[1]), f"marker mismatch {pa[:2]} vs {pb[:2]}"
        gts = [c.split(":")[0] for c in pb[9:]]  # keep phased GT only
        out.write("\t".join(pa[:9] + pa[9:] + gts) + "\n")
