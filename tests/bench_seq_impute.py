#!/usr/bin/env python3
"""Benchmark rusty-beagle (and optionally Java Beagle) on sequence-scale
imputation: one chromosome, a few thousand sequenced reference animals at ~1M
variants, and a large chip-genotyped target cohort.

Measures wall time and peak RSS per run, records output size, and deletes the
output VCF between runs so that a cohort whose full output would not fit on
disk can still be timed.

Usage:
  bench_seq_impute.py --panel DIR --sizes 250,1000,4000 [--mode impute|phase|both]
                      [--nthreads 4] [--jar beagle.jar] [--results out.tsv]

Target cohorts are generated on demand (tests/gen_seq_panel.py) and cached in
the panel directory as target_<n>.vcf.gz.
"""
import argparse, json, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "gen_seq_panel.py")


def run(cmd, log_path):
    """Run cmd, returning (wall_seconds, peak_rss_bytes, returncode).

    os.wait4 gives this child's own rusage, so the peak is per-run rather than
    the cumulative high-water mark that RUSAGE_CHILDREN would report.
    """
    exe = cmd[0] if os.path.sep in cmd[0] else shutil.which(cmd[0])
    if exe is None:
        raise FileNotFoundError(cmd[0])
    with open(log_path, "wb") as log:
        t0 = time.time()
        pid = os.posix_spawn(
            exe, cmd, os.environ,
            file_actions=[(os.POSIX_SPAWN_DUP2, log.fileno(), 1),
                          (os.POSIX_SPAWN_DUP2, log.fileno(), 2)])
        _, status, ru = os.wait4(pid, 0)
        wall = time.time() - t0
    return wall, ru.ru_maxrss * 1024, os.waitstatus_to_exitcode(status)


def ensure_target(panel, n, meta, phased, nthreads):
    """Generate target_<n>.vcf.gz (and _phased) if not already cached."""
    name = f"target_{n}"
    path = os.path.join(panel, f"{name}.vcf.gz")
    ph_path = os.path.join(panel, f"{name}_phased.vcf.gz")
    need = not os.path.exists(path) or (phased and not os.path.exists(ph_path))
    if need:
        cmd = [sys.executable, GEN, "--out-dir", panel,
               "--n-ref", str(meta["n_ref"]), "--n-targ", str(n),
               "--n-markers", str(meta["n_markers"]),
               "--chrom-len", str(meta["chrom_len"]),
               "--chip-every", str(meta["chip_every"]),
               "--chrom", str(meta["chrom"]), "--seed", str(meta["seed"]),
               "--skip-ref", "--skip-map", "--target-name", name]
        if phased:
            cmd.append("--emit-phased")
        print(f"  generating {name} ...", flush=True)
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
    return path, ph_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--sizes", required=True,
                    help="comma-separated target cohort sizes")
    ap.add_argument("--mode", default="both",
                    choices=["impute", "phase", "both"],
                    help="impute = phased target; phase = unphased target")
    ap.add_argument("--nthreads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--rusty", default=os.path.join(
        os.path.dirname(HERE), "target", "release", "rusty-beagle"))
    ap.add_argument("--jar", help="Java Beagle jar, to benchmark alongside")
    ap.add_argument("--java-xmx", default=None, help="e.g. 12g")
    ap.add_argument("--ref", default=None, help="reference file (default ref.vcf.gz)")
    ap.add_argument("--results", default=None)
    ap.add_argument("--keep-output", action="store_true")
    args = ap.parse_args()

    panel = os.path.abspath(args.panel)
    meta = json.load(open(os.path.join(panel, "panel.json")))
    ref = args.ref or os.path.join(panel, "ref.vcf.gz")
    gmap = os.path.join(panel, "plink.map")
    sizes = [int(s) for s in args.sizes.split(",")]
    modes = ["impute", "phase"] if args.mode == "both" else [args.mode]

    rows = []
    results = args.results or os.path.join(panel, "bench_results.tsv")
    if not os.path.exists(results):
        with open(results, "w") as f:
            f.write("prog\tmode\tn_targ\tn_ref\tn_markers\tn_chip\tnthreads"
                    "\twall_s\tpeak_rss_gb\tout_gb\n")

    for n in sizes:
        tgt, tgt_ph = ensure_target(panel, n, meta, "impute" in modes,
                                    args.nthreads)
        for mode in modes:
            gt = tgt_ph if mode == "impute" else tgt
            for prog in (["rusty"] + (["java"] if args.jar else [])):
                out = os.path.join(panel, f"{prog}_{mode}_{n}")
                base = [f"gt={gt}", f"ref={ref}", f"out={out}",
                        f"map={gmap}", f"nthreads={args.nthreads}"]
                if prog == "rusty":
                    cmd = [args.rusty] + base
                else:
                    java = ["java"]
                    if args.java_xmx:
                        java.append(f"-Xmx{args.java_xmx}")
                    cmd = java + ["-jar", args.jar] + base
                log = out + ".log"
                print(f"[{prog}] {mode} n_targ={n} ...", flush=True)
                wall, peak, rc = run(cmd, log)
                if rc != 0:
                    print(f"  FAILED rc={rc}; tail of {log}:", flush=True)
                    print(open(log, errors="replace").read()[-1500:], flush=True)
                    rows.append((prog, mode, n, "FAILED", rc))
                    for ext in (".vcf.gz", ".log"):
                        if os.path.exists(out + ext) and not args.keep_output:
                            os.remove(out + ext)
                    continue
                vcf = out + ".vcf.gz"
                out_gb = os.path.getsize(vcf) / 1e9 if os.path.exists(vcf) else 0.0
                print(f"  {wall:8.1f}s  peakRSS {peak/1e9:5.2f} GB  "
                      f"out {out_gb:.2f} GB", flush=True)
                with open(results, "a") as f:
                    f.write(f"{prog}\t{mode}\t{n}\t{meta['n_ref']}\t"
                            f"{meta['n_markers']}\t{meta['n_chip']}\t"
                            f"{args.nthreads}\t{wall:.1f}\t{peak/1e9:.2f}\t"
                            f"{out_gb:.2f}\n")
                if not args.keep_output:
                    for ext in (".vcf.gz",):
                        if os.path.exists(out + ext):
                            os.remove(out + ext)
    print(f"\nresults appended to {results}")


if __name__ == "__main__":
    main()
