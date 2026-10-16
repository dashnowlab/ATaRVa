#!/usr/bin/env python3
"""Compare this run's metrics (benchmark/data/latest_metrics.tsv, written by baseline_history.
append() during 04_score.py / 06_score_strts_native.py) against the checked-in accepted
baseline (benchmark/baseline_current.tsv). Exits 1 if any guarded metric dropped by more than
--tolerance, so CI can fail the job on a real accuracy regression.

Same pattern as str-truth-set-v2's own docs/truthset_release_check.md /
run_truthset_release_check.py, applied to our numbers instead of theirs: every rate here is in
[0, 1] or [0, 100] where higher is better, so one tolerance applies to all of them (a metric's
own scale is preserved between baseline_current.tsv and latest_metrics.tsv, so comparing raw
deltas is valid without normalizing).

Usage:
    python3 07_check_regression.py
    python3 07_check_regression.py --tolerance 0.02
    python3 07_check_regression.py --update-baseline   # only after confirming a drop is legitimate
"""
import argparse
import csv
import os

BENCHMARK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_BASELINE_PATH = os.path.join(BENCHMARK_DIR, "baseline_current.tsv")
DEFAULT_LATEST_PATH = os.path.join(BENCHMARK_DIR, "data", "latest_metrics.tsv")

# Percentage-point/fraction tolerance. Our n=1196 is orders of magnitude smaller than
# str-truth-set-v2's own 0.002-tolerance guards (146k+ loci there), so a rerun of the exact
# same code has more room for legitimate noise -- ATaRVa's own consensus/clustering steps
# aren't guaranteed bit-for-bit deterministic across runs/thread-scheduling.
DEFAULT_TOLERANCE = 0.01


def read_metrics_tsv(path, key_fields):
    """Read a TSV into {tuple(key_fields values): float(value)}. Skips '#'-comment lines,
    so this also reads format_baseline_tsv()'s own output back in."""
    with open(path) as f:
        lines = [line for line in f if not line.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    return {tuple(row[k] for k in key_fields): float(row["value"]) for row in rows}


def format_baseline_tsv(metrics):
    lines = [
        "# Accepted baseline for benchmark/scripts/07_check_regression.py, catalog "
        "random_2to6bp_with_extremes. A CI run's benchmark/data/latest_metrics.tsv is compared "
        "against this file; a guarded metric dropping by more than --tolerance fails the job.",
        "#",
        "# Regenerate with --update-baseline, and only after confirming a drop is a legitimate",
        "# tradeoff, not a bug -- say why in the commit message, since this is what every later",
        "# ATaRVa change is judged against.",
        "tech\tsource\tmetric\tvalue",
    ]
    for key in sorted(metrics):
        lines.append("\t".join(key) + f"\t{metrics[key]:.6f}")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default=DEFAULT_BASELINE_PATH)
    parser.add_argument("--latest", default=DEFAULT_LATEST_PATH)
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument("--update-baseline", action="store_true",
                        help="Overwrite --baseline with --latest's values instead of checking.")
    args = parser.parse_args()

    key_fields = ("tech", "source", "metric")
    latest = read_metrics_tsv(args.latest, key_fields)

    if args.update_baseline:
        with open(args.baseline, "w") as f:
            f.write(format_baseline_tsv(latest))
        print(f"Wrote {len(latest)} metrics to {args.baseline}")
        return

    baseline = read_metrics_tsv(args.baseline, key_fields)

    header = f"   {'tech':7s} {'source':18s} {'metric':40s} {'baseline':>10s} {'latest':>10s} {'delta':>9s}"
    lines = [header]
    failures = []
    for key in sorted(set(baseline) | set(latest)):
        tech, source, metric = key
        if key not in latest:
            lines.append(f"!! {tech:7s} {source:18s} {metric:40s} {baseline[key]:10.4f} {'MISSING':>10s}")
            failures.append(key)
            continue
        if key not in baseline:
            lines.append(f"   {tech:7s} {source:18s} {metric:40s} {'n/a':>10s} {latest[key]:10.4f}  (new)")
            continue

        delta = latest[key] - baseline[key]
        marker = "  "
        if delta < -args.tolerance:
            marker = "!!"
            failures.append(key)
        lines.append(f"{marker} {tech:7s} {source:18s} {metric:40s} {baseline[key]:10.4f} "
                     f"{latest[key]:10.4f} {delta:+9.4f}")

    print("\n".join(lines))

    if failures:
        print(f"\nFAIL: {len(failures)} metric(s) regressed by more than tolerance={args.tolerance} "
             f"(or are missing from this run).")
        raise SystemExit(1)

    print(f"\nOK: no metric regressed by more than tolerance={args.tolerance}.")


if __name__ == "__main__":
    main()
