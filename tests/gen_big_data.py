#!/usr/bin/env python3
"""Fast (numpy) generator of large synthetic phased VCF panels for
benchmarking rusty-beagle vs Java Beagle.  Same mosaic-of-founders model as
gen_test_data.py, vectorized."""
import argparse
import gzip
import numpy as np


def simulate(n_founders, n_haps, n_markers, region_bp, seed):
    rng = np.random.default_rng(seed)
    positions = np.sort(
        rng.choice(np.arange(10_000, region_bp), size=n_markers, replace=False))
    freq = rng.beta(0.4, 0.4, size=n_markers)
    founders = (rng.random((n_founders, n_markers)) < freq).astype(np.int8)
    # crossover mask -> segment ids -> founder per (hap, segment)
    xover = rng.random((n_haps, n_markers)) < 0.002
    seg = np.cumsum(xover, axis=1)
    max_seg = int(seg.max()) + 1
    founder_choice = rng.integers(0, n_founders, size=(n_haps, max_seg))
    haps = founders[founder_choice[np.arange(n_haps)[:, None], seg],
                    np.arange(n_markers)[None, :]]
    mut = rng.random((n_haps, n_markers)) < 0.0005
    haps = np.where(mut, 1 - haps, haps)
    return positions, haps


def write_vcf(path, chrom, positions, haps, sample_prefix, marker_idx=None):
    n_samples = haps.shape[0] // 2
    if marker_idx is None:
        marker_idx = np.arange(len(positions))
    a1 = haps[0::2, :]
    a2 = haps[1::2, :]
    with gzip.open(path, "wt", compresslevel=4) as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##source="gen_big_data.py"\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in range(n_samples):
            f.write(f"\t{sample_prefix}{s}")
        f.write("\n")
        gt_cache = {(x, y): f"\t{x}|{y}" for x in (0, 1) for y in (0, 1)}
        for m in marker_idx:
            prefix = f"{chrom}\t{positions[m]}\trs{m}\tA\tC\t.\tPASS\t.\tGT"
            row = [gt_cache[(x, y)] for x, y in zip(a1[:, m], a2[:, m])]
            f.write(prefix + "".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-ref", type=int, default=5000)
    ap.add_argument("--n-targ", type=int, default=200)
    ap.add_argument("--n-markers", type=int, default=50000)
    ap.add_argument("--region-bp", type=int, default=50_000_000)
    ap.add_argument("--keep-every", type=int, default=25)
    ap.add_argument("--n-founders", type=int, default=100)
    ap.add_argument("--seed", type=int, default=99)
    ap.add_argument("--chrom", default="20")
    args = ap.parse_args()

    n_haps = 2 * (args.n_ref + args.n_targ)
    positions, haps = simulate(args.n_founders, n_haps, args.n_markers,
                               args.region_bp, args.seed)
    keep = np.arange(0, args.n_markers, args.keep_every)
    write_vcf(f"{args.out_dir}/ref.vcf.gz", args.chrom, positions,
              haps[: 2 * args.n_ref], "REF")
    write_vcf(f"{args.out_dir}/target.vcf.gz", args.chrom, positions,
              haps[2 * args.n_ref:], "TARG", marker_idx=keep)
    print(f"wrote {args.n_ref} ref samples x {args.n_markers} markers; "
          f"{args.n_targ} targ samples x {len(keep)} markers")


if __name__ == "__main__":
    main()
