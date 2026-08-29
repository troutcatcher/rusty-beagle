"""Generate a livestock-shaped panel: many animals, one chromosome, few SNPs.

Usage: gen_large_cohort.py <n_ref> <n_targ> <n_markers> <chrom_len_bp> <outdir>
e.g.   gen_large_cohort.py 200000 5000 1300 158000000 /tmp/big
then   java -jar bref3.jar ref.vcf > ref.bref3


Reference haplotypes are mosaics of a small founder pool, so the panel has
realistic LD (which is what makes imputation and the sequence coder behave
like they would on real data).
"""
import numpy as np, sys, time

n_ref     = int(sys.argv[1])          # reference animals
n_targ    = int(sys.argv[2])          # target animals
n_markers = int(sys.argv[3])          # SNPs on the chromosome
chrom_len = int(sys.argv[4])          # chromosome length in bp
outdir    = sys.argv[5]
seed      = 20250829

rng = np.random.default_rng(seed)
N_FOUNDERS = 60
# ~1 recombination per 100 Mb per haplotype, expressed per marker interval
RECOMB_PER_MARKER = (chrom_len / n_markers) / 1e8

pos = np.sort(rng.choice(np.arange(1000, chrom_len, 50), size=n_markers, replace=False))
maf = rng.uniform(0.05, 0.5, size=n_markers)
founders = (rng.random((N_FOUNDERS, n_markers)) < maf).astype(np.int8)

def mosaic(n_haps, chunk=8000):
    """(n_markers, n_haps) int8 allele matrix, marker-major like a VCF."""
    out = np.empty((n_markers, n_haps), dtype=np.int8)
    for lo in range(0, n_haps, chunk):
        hi = min(lo + chunk, n_haps)
        h = hi - lo
        switch = rng.random((h, n_markers)) < RECOMB_PER_MARKER
        switch[:, 0] = True
        idx = np.where(switch, np.arange(n_markers, dtype=np.int32), 0)
        np.maximum.accumulate(idx, axis=1, out=idx)
        who = rng.integers(0, N_FOUNDERS, size=(h, n_markers), dtype=np.int8)
        who = np.take_along_axis(who, idx, axis=1)
        out[:, lo:hi] = founders[who, np.arange(n_markers)].T
    return out

# "a|b\t" for each of the 4 phased genotypes, so a line body is one memcpy
GT = np.frombuffer(b"0|0\t0|1\t1|0\t1|1\t", dtype=np.uint8).reshape(4, 4)

def write_vcf(path, prefix, alleles, markers, phased=True):
    n_s = alleles.shape[1] // 2
    ids = "\t".join(f"{prefix}{i}" for i in range(n_s))
    with open(path, "wb") as f:
        f.write(b"##fileformat=VCFv4.2\n##contig=<ID=1>\n")
        f.write(("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + ids + "\n").encode())
        sep = b"|" if phased else b"/"
        for j, m in enumerate(markers):
            row = alleles[m]
            code = (row[0::2].astype(np.uint8) << 1) | row[1::2].astype(np.uint8)
            body = GT[code].tobytes()[:-1]          # drop trailing tab
            if not phased:
                body = body.replace(b"|", sep)
            f.write(b"1\t%d\trs%d\tA\tC\t.\tPASS\t.\tGT\t" % (pos[m], m) + body + b"\n")

t0 = time.time()
ref = mosaic(2 * n_ref)
print(f"ref mosaic {ref.shape} in {time.time()-t0:.1f}s", flush=True)
write_vcf(f"{outdir}/ref.vcf", "REF", ref, range(n_markers))
print(f"ref.vcf written in {time.time()-t0:.1f}s", flush=True)
del ref

targ = mosaic(2 * n_targ)
# target animals are genotyped on a low-density subset of the same chip
chip = list(range(0, n_markers, 10))
write_vcf(f"{outdir}/target.vcf", "TARG", targ, chip)
write_vcf(f"{outdir}/target_unph.vcf", "TARG", targ, chip, phased=False)
print(f"target: {len(chip)} of {n_markers} markers, {n_targ} animals", flush=True)
print(f"total {time.time()-t0:.1f}s", flush=True)
