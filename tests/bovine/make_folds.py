"""Builds sire-family-aware cross-validation folds for the synbreedData
cattle panel and writes the per-fold VCFs.

Whole paternal half-sib families are assigned to one fold, so a target
animal never has a half-sib in its reference panel -- random splits leak
family haplotypes and overstate imputation accuracy.

Per fold k:
  ref_k.vcf.gz    all markers, the other ~400 bulls, unphased (missing kept)
  targ_k.vcf.gz   every Nth marker per chromosome, fold bulls, unphased
  truth_k.npz     masked-marker genotypes of the fold bulls + marker metadata
"""
import argparse, gzip, sys
import numpy as np

def read_geno(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        snps = header[1:]
        ids, rows = [], []
        for line in f:
            p = line.rstrip("\n").split("\t")
            ids.append(p[0])
            rows.append([-1 if v == "NA" else int(v) for v in p[1:]])
    return ids, snps, np.array(rows, dtype=np.int8)

def read_map(path):
    out = []
    with open(path) as f:
        f.readline()
        for line in f:
            snp, chrom, pos_mb = line.split("\t")
            out.append((snp, int(chrom), float(pos_mb)))
    return out

def clean_map(snps, mp):
    """Sort by (chr, bp); drop markers whose rounded bp collides."""
    idx = {s: i for i, s in enumerate(snps)}
    rows = [(c, max(1, round(p * 1e6)), s, idx[s]) for (s, c, p) in mp]
    rows.sort()
    keep, last = [], None
    for c, bp, s, i in rows:
        if (c, bp) == last:
            continue
        last = (c, bp)
        keep.append((c, bp, s, i))
    return keep  # list of (chr, bp, snp, column-index) in output order

def sire_folds(ids, ped_path, n_folds, seed):
    sire = {}
    with open(ped_path) as f:
        f.readline()
        for line in f:
            p = line.split("\t")
            sire[p[0]] = p[1]
    fams = {}
    for a in ids:
        key = sire.get(a, "0")
        if key == "0":
            key = "self_" + a  # unknown sire: own singleton family
        fams.setdefault(key, []).append(a)
    rng = np.random.default_rng(seed)
    order = sorted(fams)
    rng.shuffle(order)
    folds = [[] for _ in range(n_folds)]
    for fam in order:  # greedy: put each family in the smallest fold
        k = min(range(n_folds), key=lambda j: len(folds[j]))
        folds[k].extend(fams[fam])
    return folds

GT = {-1: "./.", 0: "0/0", 1: "0/1", 2: "1/1"}

def write_vcf(path, sample_ids, geno, markers):
    """markers: (chr, bp, snp, col) rows; geno: full matrix, rows = samples."""
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in range(1, 30):
            f.write(f"##contig=<ID={c}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(sample_ids) + "\n")
        for c, bp, snp, col in markers:
            gts = "\t".join(GT[int(v)] for v in geno[:, col])
            f.write(f"{c}\t{bp}\t{snp}\tA\tC\t.\tPASS\t.\tGT\t{gts}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--keep-every", type=int, default=5,
                    help="LD panel keeps every Nth marker per chromosome")
    ap.add_argument("--seed", type=int, default=17)
    a = ap.parse_args()

    ids, snps, geno = read_geno(f"{a.data_dir}/geno.tsv")
    markers = clean_map(snps, read_map(f"{a.data_dir}/map.tsv"))
    print(f"{len(ids)} animals, {len(markers)} of {len(snps)} markers after map cleaning")

    ld, cur_chr, i_on_chr = [], None, 0
    for row in markers:
        if row[0] != cur_chr:
            cur_chr, i_on_chr = row[0], 0
        if i_on_chr % a.keep_every == 0:
            ld.append(row)
        i_on_chr += 1
    ld_set = {r[2] for r in ld}
    masked = [r for r in markers if r[2] not in ld_set]
    print(f"LD panel {len(ld)} markers; scoring on {len(masked)} masked markers")

    folds = sire_folds(ids, f"{a.data_dir}/ped.tsv", a.folds, a.seed)
    id_row = {s: i for i, s in enumerate(ids)}
    for k, fold in enumerate(folds):
        targ = sorted(fold)
        ref = sorted(set(ids) - set(fold))
        tr = np.array([id_row[s] for s in targ])
        rr = np.array([id_row[s] for s in ref])
        write_vcf(f"{a.out_dir}/ref_{k}.vcf.gz", ref, geno[rr], markers)
        write_vcf(f"{a.out_dir}/targ_{k}.vcf.gz", targ, geno[tr], ld)
        truth = geno[np.ix_(tr, [m[3] for m in masked])]
        # reference-panel allele frequency of the ALT (B) allele per masked marker
        sub = geno[np.ix_(rr, [m[3] for m in masked])].astype(np.float64)
        sub[sub < 0] = np.nan
        af = np.nanmean(sub, axis=0) / 2.0
        np.savez_compressed(
            f"{a.out_dir}/truth_{k}.npz",
            truth=truth,
            samples=np.array(targ),
            chrom=np.array([m[0] for m in masked]),
            pos=np.array([m[1] for m in masked]),
            snp=np.array([m[2] for m in masked]),
            ref_af=af,
        )
        print(f"fold {k}: {len(targ)} targets / {len(ref)} reference")

if __name__ == "__main__":
    main()
