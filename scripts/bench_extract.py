#!/usr/bin/env python3
"""Benchmark cfextract.extract_features() for wall time and peak memory.

Reports aggregate statistics across all samples in the manifest.

Example:
  python scripts/bench_extract.py
  python scripts/bench_extract.py --repeat 3
"""

import argparse
import csv
import json
import statistics
import subprocess
import sys
import textwrap
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", default="data/samples.csv")
    p.add_argument("--contig", default="chr21")
    p.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="extractions per sample; stats pooled across all samples (default: 1)",
    )
    return p.parse_args()


def load_bam_paths(path: Path) -> list[str]:
    base = path.parent
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return [str((base / row["bam_path"]).resolve()) for row in rows]


# Worker runs in a fresh subprocess per measurement so the OS allocator
# cannot recycle pages from a prior call, giving an accurate RSS delta.
_WORKER = textwrap.dedent("""\
    import json, sys, time, cfextract, psutil
    bam, contig = sys.argv[1], sys.argv[2]
    proc = psutil.Process()
    rss0 = proc.memory_info().rss
    t0 = time.perf_counter()
    cfextract.extract_features(bam, contig)
    elapsed = time.perf_counter() - t0
    rss1 = proc.memory_info().rss
    print(json.dumps({"time_s": elapsed, "mem_mb": (rss1 - rss0) / 1024**2}))
""")
def benchmark(bam_abs: str, contig: str) -> dict[str, float]:
    result = subprocess.run(
        [sys.executable, "-c", _WORKER, bam_abs, contig],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout.strip())


def main() -> None:
    args = parse_args()
    bam_paths = load_bam_paths(Path(args.manifest))
    n = len(bam_paths)

    all_times: list[float] = []
    all_mems: list[float] = []

    for i, bam in enumerate(bam_paths, 1):
        print(f"[{i}/{n}] {Path(bam).name} ...", end=" ", flush=True, file=sys.stderr)
        for _ in range(args.repeat):
            r = benchmark(bam, args.contig)
            all_times.append(r["time_s"])
            all_mems.append(r["mem_mb"])
        print("done", file=sys.stderr)

    t_mean = statistics.mean(all_times)
    t_med = statistics.median(all_times)
    t_std = statistics.stdev(all_times) if len(all_times) > 1 else 0.0
    m_mean = statistics.mean(all_mems)
    m_med = statistics.median(all_mems)
    m_std = statistics.stdev(all_mems) if len(all_mems) > 1 else 0.0

    print(
        f"cfextract.extract_features()"
        f" — contig={args.contig}  n={n}  repeat={args.repeat}"
    )
    print(f"  time:   {t_mean:.2f} ± {t_std:.2f} s   (median {t_med:.2f} s)")
    print(f"  memory: {m_mean:.1f} ± {m_std:.1f} MB  (median {m_med:.1f} MB)")


if __name__ == "__main__":
    main()
