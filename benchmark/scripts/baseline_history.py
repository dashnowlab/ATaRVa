#!/usr/bin/env python3
"""Shared helper: append a scoring run's metrics to benchmark/baseline_history.tsv, and record
them for this run in benchmark/data/latest_metrics.tsv (07_check_regression.py's input).

baseline_history.tsv is one flat, append-only TSV (date, atarva_version, tech, catalog, source,
metric, value) rather than a single per-tech snapshot file, so results from every ATaRVa
version ever benchmarked stay around for comparison -- the point is tracking history across
releases, not just "the current number". Load it with pandas / a pivot table / `awk` to
compare, e.g.:

    awk -F'\t' '$5=="BASEPAIR/f1"' benchmark/baseline_history.tsv

latest_metrics.tsv holds only this run's rows (tech, source, metric, value -- no date/version/
catalog, since it's scoped to a single run) and is appended to across the handful of scoring
calls in one CI/local run (04_score.py x2 techs + 06_score_strts_native.py x2 techs);
07_check_regression.py reads it back to compare against baseline_current.tsv. Delete it
(`rm -f benchmark/data/latest_metrics.tsv`) before starting a fresh scoring pass, same as any
other file under data/ -- it isn't reset automatically, so a stale file from a previous/partial
run would otherwise leak into a new regression check.
"""
import csv
import datetime
import os

BENCHMARK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HISTORY_PATH = os.path.join(BENCHMARK_DIR, "baseline_history.tsv")
LATEST_METRICS_PATH = os.path.join(BENCHMARK_DIR, "data", "latest_metrics.tsv")

COLUMNS = ["date", "atarva_version", "tech", "catalog", "source", "metric", "value"]
LATEST_COLUMNS = ["tech", "source", "metric", "value"]


def append(atarva_version, tech, catalog, source, metrics, history_path=HISTORY_PATH,
          latest_metrics_path=LATEST_METRICS_PATH):
    """Append one row per metric to the history log and to this run's latest-metrics file.

    Args:
        atarva_version (str): e.g. "0.7.1+ext0.01" -- from `atarva[-ext] --version`.
        tech (str): "pacbio" or "ont".
        catalog (str): which locus catalog this run used, e.g. "random_2to6bp_with_extremes".
        source (str): which comparison produced these metrics, e.g. "aardvark_basepair" or
            "strts_native".
        metrics (dict): metric name -> value.
        history_path (str): defaults to benchmark/baseline_history.tsv.
        latest_metrics_path (str): defaults to benchmark/data/latest_metrics.tsv.
    """
    file_exists = os.path.exists(history_path)
    date = datetime.date.today().isoformat()
    with open(history_path, "a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        if not file_exists:
            writer.writerow(COLUMNS)
        for metric in sorted(metrics):
            writer.writerow([date, atarva_version, tech, catalog, source, metric, f"{metrics[metric]:.6f}"])
    print(f"Appended {len(metrics)} rows to {history_path} "
         f"(version={atarva_version}, tech={tech}, catalog={catalog}, source={source})")

    latest_file_exists = os.path.exists(latest_metrics_path)
    with open(latest_metrics_path, "a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        if not latest_file_exists:
            writer.writerow(LATEST_COLUMNS)
        for metric in sorted(metrics):
            writer.writerow([tech, source, metric, f"{metrics[metric]:.6f}"])
