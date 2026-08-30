#!/usr/bin/env python3
"""Generate a sequence-imputation panel: one chromosome, ~1M sequence variants,
a few thousand sequenced reference animals, and a large chip-genotyped target
cohort.

This is the "impute to sequence" shape used in livestock genomics: a small
whole-genome-sequenced reference (thousands of animals, millions of variants,
a site-frequency spectrum dominated by rare alleles) and a very large target
cohort carrying only the chip subset.

Usage:
  gen_seq_panel.py --out-dir DIR --n-ref 3500 --n-targ 50000 \
                   --n-markers 1000000 --chrom-len 158000000 --chip-every 40

Writes (BGZF-compressed, as real pipelines do):
  ref.vcf.gz          n_ref phased animals x n_markers sequence variants
  target.vcf.gz       n_targ unphased animals x the chip subset
  target_phased.vcf.gz  same, phased        (--emit-phased)
  plink.map           PLINK genetic map for the chromosome
  panel.json          marker/sample counts, for the benchmark harness

Reference haplotypes are mosaics of a founder pool, so the panel carries
realistic LD -- which is what makes imputation, the sequence coder and the
rare-variant (second-stage) phasing behave like they would on real data.
"""
import argparse, json, os, struct, sys, time, zlib
import multiprocessing as mp
import numpy as np

# ---------------------------------------------------------------- BGZF output
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
MAX_BLOCK = 65000          # payload per BGZF block (limit is 65535 compressed)


def _deflate(payload, level):
    co = zlib.compressobj(level, zlib.DEFLATED, -15)
    data = co.compress(payload) + co.flush()
    bsize = len(data) + 25 + 1
    if bsize > 65536:                      # incompressible: fall back to store
        co = zlib.compressobj(0, zlib.DEFLATED, -15)
        data = co.compress(payload) + co.flush()
        bsize = len(data) + 25 + 1
    return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
            + struct.pack("<H", bsize - 1) + data
            + struct.pack("<II", zlib.crc32(payload), len(payload)))


class BgzfWriter:
    """Buffered BGZF writer; blocks are compressed on a worker pool."""

    def __init__(self, path, level=4, procs=None):
        self.f = open(path, "wb")
        self.level = level
        self.buf = bytearray()
        self.pending = []
        self.pool = mp.Pool(procs or max(1, (os.cpu_count() or 2)))
        self.batch = 256                    # blocks dispatched per pool round

    def write(self, chunk):
        self.buf += chunk
        if len(self.buf) >= MAX_BLOCK * self.batch:
            self._flush_blocks(keep_tail=True)

    def _flush_blocks(self, keep_tail):
        n = len(self.buf) // MAX_BLOCK if keep_tail else -(-len(self.buf) // MAX_BLOCK)
        if n <= 0:
            return
        blocks = [bytes(self.buf[i * MAX_BLOCK:(i + 1) * MAX_BLOCK]) for i in range(n)]
        del self.buf[:n * MAX_BLOCK]
        for out in self.pool.starmap(_deflate, [(b, self.level) for b in blocks]):
            self.f.write(out)

    def close(self):
        self._flush_blocks(keep_tail=False)
        self.f.write(BGZF_EOF)
        self.f.close()
        self.pool.close()
        self.pool.join()


# ------------------------------------------------------------------ simulation
def site_frequency_spectrum(rng, n, f_min, f_max):
    """Neutral-ish SFS: density proportional to 1/f over [f_min, f_max]."""
    u = rng.random(n)
    return f_min * (f_max / f_min) ** u


def marker_positions(rng, n_markers, chrom_len):
    """Strictly increasing positions with mean spacing chrom_len/n_markers."""
    mean_gap = chrom_len / n_markers
    gaps = rng.integers(1, max(2, int(2 * mean_gap)), size=n_markers, dtype=np.int64)
    return 1000 + np.cumsum(gaps)


# "a|b\t" for each of the 4 phased genotypes, so a line body is one memcpy
GT_PH = np.frombuffer(b"0|0\t0|1\t1|0\t1|1\t", dtype=np.uint8).reshape(4, 4)
GT_UN = np.frombuffer(b"0/0\t0/1\t0/1\t1/1\t", dtype=np.uint8).reshape(4, 4)


def vcf_header(chrom, chrom_len, prefix, n_samples, offset=0):
    ids = "\t".join(f"{prefix}{offset + i}" for i in range(n_samples))
    return ("##fileformat=VCFv4.2\n"
            f"##contig=<ID={chrom},length={chrom_len}>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + ids + "\n").encode()


def stream_panel(path, chrom, chrom_len, pos, founders, markers, n_haps, recomb,
                 rng, prefix, phased, level, block_haps=50_000_000, offset=0):
    """Write a mosaic-of-founders VCF, generating markers in blocks so that
    peak memory stays bounded regardless of cohort size."""
    n_m = len(markers)
    B = max(1, min(n_m, block_haps // max(1, n_haps)))
    w = BgzfWriter(path, level=level)
    w.write(vcf_header(chrom, chrom_len, prefix, n_haps // 2, offset))
    table = GT_PH if phased else GT_UN
    n_f = founders.shape[0]
    who = rng.integers(0, n_f, size=n_haps, dtype=np.int16)   # carried across blocks
    ar = np.arange(B, dtype=np.int32)
    t0 = time.time()
    for lo in range(0, n_m, B):
        hi = min(lo + B, n_m)
        b = hi - lo
        sw = rng.random((n_haps, b)) < recomb
        idx = np.where(sw, ar[:b], -1)
        np.maximum.accumulate(idx, axis=1, out=idx)
        pick = rng.integers(0, n_f, size=(n_haps, b), dtype=np.int16)
        pick = np.take_along_axis(pick, np.maximum(idx, 0), axis=1)
        chosen = np.where(idx >= 0, pick, who[:, None])
        who = chosen[:, -1].copy()
        cols = markers[lo:hi]
        alleles = founders[chosen, cols[None, :]]             # (n_haps, b)
        del sw, idx, pick, chosen
        codes = ((alleles[0::2, :].astype(np.uint8) << 1)
                 | alleles[1::2, :].astype(np.uint8))         # (n_samples, b)
        del alleles
        out = []
        for j in range(b):
            m = cols[j]
            body = table[codes[:, j]].tobytes()[:-1]          # drop trailing tab
            out.append(b"%s\t%d\trs%d\tA\tC\t.\tPASS\t.\tGT\t"
                       % (chrom.encode(), pos[m], m) + body + b"\n")
        w.write(b"".join(out))
        del codes, out
        done = hi
        print(f"    {os.path.basename(path)}: {done}/{n_m} markers "
              f"({time.time()-t0:.0f}s)", flush=True)
    w.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-ref", type=int, default=3500)
    ap.add_argument("--n-targ", type=int, default=50000)
    ap.add_argument("--n-markers", type=int, default=1_000_000)
    ap.add_argument("--chrom-len", type=int, default=158_000_000)
    ap.add_argument("--chip-every", type=int, default=40,
                    help="target animals carry every Nth sequence variant")
    ap.add_argument("--n-founders", type=int, default=200)
    ap.add_argument("--cm-per-mb", type=float, default=1.0)
    ap.add_argument("--chrom", default="1")
    ap.add_argument("--seed", type=int, default=20260830)
    ap.add_argument("--level", type=int, default=4, help="BGZF deflate level")
    ap.add_argument("--emit-phased", action="store_true",
                    help="also write a phased target (imputation-only runs)")
    ap.add_argument("--skip-ref", action="store_true")
    ap.add_argument("--skip-target", action="store_true")
    ap.add_argument("--skip-map", action="store_true")
    ap.add_argument("--target-name", default="target",
                    help="basename for the target VCF(s)")
    ap.add_argument("--batch-size", type=int, default=0,
                    help="split the target cohort into disjoint batch files "
                         "of this many animals (0 = one file)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    n_m = args.n_markers

    t0 = time.time()
    pos = marker_positions(rng, n_m, args.chrom_len)
    chrom_len = int(pos[-1]) + 1000
    freq = site_frequency_spectrum(rng, n_m, 1.0 / args.n_founders, 0.5)
    founders = (rng.random((args.n_founders, n_m)) < freq).astype(np.int8)
    print(f"founders {founders.shape} in {time.time()-t0:.0f}s", flush=True)

    # PLINK genetic map for the chromosome
    cm = pos * (args.cm_per_mb / 1e6)
    if not args.skip_map:
        with open(f"{args.out_dir}/plink.map", "w") as f:
            f.write("".join(f"{args.chrom}\trs{m}\t{cm[m]:.8f}\t{pos[m]}\n"
                            for m in range(n_m)))
        print(f"plink.map written ({cm[-1]:.1f} cM)", flush=True)

    all_markers = np.arange(n_m, dtype=np.int64)
    chip = np.arange(0, n_m, args.chip_every, dtype=np.int64)
    # 1 crossover per 100 Mb per haplotype, per marker interval
    recomb_ref = (args.chrom_len / n_m) / 1e8
    recomb_chip = recomb_ref * args.chip_every

    if not args.skip_ref:
        stream_panel(f"{args.out_dir}/ref.vcf.gz", args.chrom, chrom_len, pos,
                     founders, all_markers, 2 * args.n_ref, recomb_ref, rng,
                     "REF", True, args.level)
        print(f"ref done {time.time()-t0:.0f}s", flush=True)

    if not args.skip_target:
        st = np.random.default_rng(args.seed + 1)
        bs = args.batch_size or args.n_targ
        n_b = -(-args.n_targ // bs)
        for b in range(n_b):
            k = min(bs, args.n_targ - b * bs)
            tag = args.target_name if n_b == 1 else f"{args.target_name}_b{b:02d}"
            stream_panel(f"{args.out_dir}/{tag}.vcf.gz", args.chrom, chrom_len,
                         pos, founders, chip, 2 * k, recomb_chip, st, "TARG",
                         False, args.level, offset=b * bs)
            if args.emit_phased:
                st2 = np.random.default_rng(args.seed + 1000 + b)
                stream_panel(f"{args.out_dir}/{tag}_phased.vcf.gz", args.chrom,
                             chrom_len, pos, founders, chip, 2 * k, recomb_chip,
                             st2, "TARG", True, args.level, offset=b * bs)

    meta = "panel.json" if args.skip_target else f"{args.target_name}.json"
    with open(f"{args.out_dir}/{meta}", "w") as f:
        json.dump({"n_ref": args.n_ref, "n_targ": args.n_targ,
                   "n_markers": n_m, "n_chip": int(len(chip)),
                   "chip_every": args.chip_every, "chrom": args.chrom,
                   "chrom_len": chrom_len, "cM": float(cm[-1]),
                   "seed": args.seed}, f, indent=2)
    print(f"total {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
