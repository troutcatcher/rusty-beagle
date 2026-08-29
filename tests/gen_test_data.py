#!/usr/bin/env python3
"""Generates synthetic phased VCF test data for validating rusty-beagle
against Java Beagle.

Simulates haplotypes as mosaics of founder haplotypes (crossovers +
mutation), giving realistic LD.  Produces:
  ref.vcf.gz    - phased reference panel, all markers
  target.vcf.gz - phased target samples, every k-th marker retained
  truth.vcf.gz  - phased target samples, all markers (for accuracy checks)
"""
import argparse
import gzip
import random


def simulate(n_founders, n_haps, n_markers, region_bp, seed, multi_prop=0.0):
    rng = random.Random(seed)
    positions = sorted(rng.sample(range(10_000, region_bp), n_markers))
    # allele frequency per marker (founder haplotype alleles)
    n_alleles = []
    for m in range(n_markers):
        if multi_prop > 0 and rng.random() < multi_prop:
            n_alleles.append(rng.choice([3, 4]))
        else:
            n_alleles.append(2)
    founders = []
    for _ in range(n_founders):
        hap = []
        for m in range(n_markers):
            f = rng.betavariate(0.4, 0.4)
            hap.append(0 if rng.random() > f else rng.randrange(n_alleles[m]))
        founders.append(hap)
    # correlate founders: copy blocks between founders
    for f in range(1, n_founders):
        src = rng.randrange(f)
        start = rng.randrange(n_markers)
        end = min(n_markers, start + rng.randrange(n_markers // 4 + 1))
        founders[f][start:end] = founders[src][start:end]

    haps = []
    for _ in range(n_haps):
        hap = []
        cur = rng.randrange(n_founders)
        for m in range(n_markers):
            if rng.random() < 0.002:  # crossover
                cur = rng.randrange(n_founders)
            allele = founders[cur][m]
            if rng.random() < 0.0005:  # mutation
                allele = rng.randrange(n_alleles[m])
            hap.append(allele)
        haps.append(hap)
    return positions, n_alleles, haps


ALTS = ["C", "G", "T"]


def write_vcf(path, chrom, positions, n_alleles, haps, sample_prefix,
              marker_idx=None, haploid_samples=(), with_end_info=False):
    n_samples = len(haps) // 2
    if marker_idx is None:
        marker_idx = range(len(positions))
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##source="gen_test_data.py"\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in range(n_samples):
            f.write(f"\t{sample_prefix}{s}")
        f.write("\n")
        for m in marker_idx:
            alt = ",".join(ALTS[: n_alleles[m] - 1])
            marker_id = f"rs{m}" if m % 3 else "."
            info = f"END={positions[m]}" if (with_end_info and m % 17 == 0) else "."
            f.write(f"{chrom}\t{positions[m]}\t{marker_id}\tA\t{alt}\t.\tPASS\t{info}\tGT")
            for s in range(n_samples):
                a1 = haps[2 * s][m]
                a2 = haps[2 * s + 1][m]
                if s in haploid_samples:
                    f.write(f"\t{a1}")
                else:
                    f.write(f"\t{a1}|{a2}")
            f.write("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-ref", type=int, default=190)
    ap.add_argument("--n-targ", type=int, default=10)
    ap.add_argument("--n-markers", type=int, default=1000)
    ap.add_argument("--region-bp", type=int, default=2_000_000)
    ap.add_argument("--keep-every", type=int, default=10)
    ap.add_argument("--n-founders", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--chrom", default="22")
    ap.add_argument("--multi-prop", type=float, default=0.0,
                    help="proportion of multi-allelic markers")
    ap.add_argument("--haploid-targ", type=int, default=0,
                    help="number of haploid target samples")
    ap.add_argument("--with-end-info", action="store_true")
    args = ap.parse_args()

    n_haps = 2 * (args.n_ref + args.n_targ)
    positions, n_alleles, haps = simulate(
        args.n_founders, n_haps, args.n_markers, args.region_bp, args.seed,
        args.multi_prop)
    ref_haps = haps[: 2 * args.n_ref]
    targ_haps = haps[2 * args.n_ref:]

    keep = [m for m in range(args.n_markers) if m % args.keep_every == 0]
    haploid = set(range(args.haploid_targ))
    write_vcf(f"{args.out_dir}/ref.vcf.gz", args.chrom, positions, n_alleles,
              ref_haps, "REF", with_end_info=args.with_end_info)
    write_vcf(f"{args.out_dir}/target.vcf.gz", args.chrom, positions, n_alleles,
              targ_haps, "TARG", marker_idx=keep, haploid_samples=haploid)
    write_vcf(f"{args.out_dir}/truth.vcf.gz", args.chrom, positions, n_alleles,
              targ_haps, "TARG", haploid_samples=haploid)
    print(f"wrote {args.out_dir}/ref.vcf.gz ({args.n_ref} samples, "
          f"{args.n_markers} markers), target.vcf.gz ({args.n_targ} samples, "
          f"{len(keep)} markers)")


if __name__ == "__main__":
    main()
