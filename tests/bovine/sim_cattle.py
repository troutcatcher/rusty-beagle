"""Coalescent simulation of a cattle-like population for imputation-accuracy
experiments, using the Holstein effective-population-size trajectory of
MacLeod et al. (2013, Mol Biol Evol 30:2209), as encoded by AlphaSimR's
runMacs(species="CATTLE"): Ne=90 at present, growing to 62,000 at 6,000
generations ago; mu=2.5e-8, rec=1e-8 (~1 cM/Mb, the usual cattle average).

Outputs a phased reference VCF per requested panel size (nested subsets),
an unphased low-density target VCF, and a truth npz for the masked markers.
Chip ascertainment: HD sites are drawn evenly spaced among MAF>=0.05
variants (mimicking array design); the LD panel keeps every Nth HD site.
"""
import argparse, gzip
import numpy as np
import msprime

HIST_GEN = [3, 6, 12, 18, 24, 30, 36, 48, 124, 300, 1000, 6000]
HIST_NE = [120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000]

def cattle_demography():
    d = msprime.Demography()
    d.add_population(name="cattle", initial_size=90)
    for g, ne in zip(HIST_GEN, HIST_NE):
        d.add_population_parameters_change(time=g, initial_size=ne)
    return d

def write_ref_vcf(path, gm, sample_lo, sample_hi, chrom, positions):
    n = sample_hi - sample_lo
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=%d>\n" % chrom)
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(f"R{j}" for j in range(n)) + "\n")
        for i, pos in enumerate(positions):
            row = gm[i, 2 * sample_lo:2 * sample_hi]
            gts = "\t".join(f"{row[2*j]}|{row[2*j+1]}" for j in range(n))
            f.write(f"{chrom}\t{pos}\trs{i}\tA\tC\t.\tPASS\t.\tGT\t{gts}\n")

def write_targ_vcf(path, gm, sample_lo, sample_hi, chrom, positions, site_idx):
    n = sample_hi - sample_lo
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=%d>\n" % chrom)
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(f"T{j}" for j in range(n)) + "\n")
        for i in site_idx:
            row = gm[i, 2 * sample_lo:2 * sample_hi]
            gts = "\t".join(
                "/".join(sorted(str(a) for a in (row[2*j], row[2*j+1])))
                for j in range(n))
            f.write(f"{chrom}\t{positions[i]}\trs{i}\tA\tC\t.\tPASS\t.\tGT\t{gts}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--length-mb", type=float, default=100.0)
    ap.add_argument("--hd-per-mb", type=float, default=8.2,
                    help="HD marker density (50k-chip density on the 2.7Gb genome)")
    ap.add_argument("--keep-every", type=int, default=10)
    ap.add_argument("--ref-sizes", type=int, nargs="+", default=[400, 5000, 20000])
    ap.add_argument("--n-targ", type=int, default=250)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    L = int(a.length_mb * 1e6)
    n_dip = max(a.ref_sizes) + a.n_targ
    print(f"simulating {n_dip} diploids over {a.length_mb} Mb ...", flush=True)
    ts = msprime.sim_ancestry(samples=n_dip, demography=cattle_demography(),
                              sequence_length=L, recombination_rate=1e-8,
                              random_seed=a.seed)
    ts = msprime.sim_mutations(ts, rate=2.5e-8, random_seed=a.seed + 1,
                               model=msprime.BinaryMutationModel())
    print(f"{ts.num_sites} segregating sites", flush=True)

    # Ascertain the HD chip: evenly spaced among common variants. Allele
    # counts come from tree topology (num_samples below the mutation node),
    # which avoids decoding genotypes at every segregating site; only the
    # chosen chip sites are ever decoded.
    n_hd = int(a.length_mb * a.hd_per_mb)
    n_haps = 2 * n_dip
    common_pos, common_site = [], []
    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) != 1:
                continue  # recurrent mutation: skip
            ac = tree.num_samples(site.mutations[0].node)
            if 0.05 <= ac / n_haps <= 0.95:
                common_pos.append(site.position)
                common_site.append(site.id)
    common_pos = np.array(common_pos)
    grid = np.linspace(common_pos[0], common_pos[-1], n_hd)
    chosen = np.unique(np.searchsorted(common_pos, grid).clip(0, len(common_pos) - 1))
    print(f"{len(common_pos)} common variants; HD chip {len(chosen)} sites", flush=True)

    import tskit
    var = tskit.Variant(ts)
    gm_rows, positions = [], []
    for ci in chosen:
        var.decode(int(common_site[ci]))
        gm_rows.append(var.genotypes.astype(np.int8))
        positions.append(int(common_pos[ci]) + 1)
    gm = np.array(gm_rows)
    # collapse duplicate integer positions from rounding
    keep = [0]
    for i in range(1, len(positions)):
        if positions[i] != positions[keep[-1]]:
            keep.append(i)
    gm, positions = gm[keep], [positions[i] for i in keep]
    print(f"{len(positions)} HD markers after position rounding", flush=True)

    ld_idx = list(range(0, len(positions), a.keep_every))
    masked_idx = [i for i in range(len(positions)) if i not in set(ld_idx)]

    n_targ = a.n_targ
    targ_lo, targ_hi = 0, n_targ  # targets first; refs are nested after them
    for rs in a.ref_sizes:
        write_ref_vcf(f"{a.out_dir}/simref_{rs}.vcf.gz", gm,
                      n_targ, n_targ + rs, 1, positions)
    write_targ_vcf(f"{a.out_dir}/simtarg.vcf.gz", gm, targ_lo, targ_hi,
                   1, positions, ld_idx)

    truth = np.zeros((n_targ, len(masked_idx)), dtype=np.int8)
    for col, i in enumerate(masked_idx):
        row = gm[i, 2 * targ_lo:2 * targ_hi]
        truth[:, col] = row[0::2] + row[1::2]
    ref_af = np.array([gm[i, 2 * n_targ:].mean() for i in masked_idx])
    np.savez_compressed(
        f"{a.out_dir}/simtruth.npz", truth=truth,
        samples=np.array([f"T{j}" for j in range(n_targ)]),
        chrom=np.array([1] * len(masked_idx)),
        pos=np.array([positions[i] for i in masked_idx]),
        snp=np.array([f"rs{i}" for i in masked_idx]),
        ref_af=ref_af)
    print("done", flush=True)

if __name__ == "__main__":
    main()
