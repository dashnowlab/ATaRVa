# Reused from broadinstitute/str-truth-set

These files are copied, unmodified, from the [broadinstitute/str-truth-set](https://github.com/broadinstitute/str-truth-set)
repository (MIT licensed), commit `main` as of 2026-09-11:

- `tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py`
- `tool_comparison/scripts/add_tool_results_columns.py`
- `tool_comparison/scripts/add_concordance_columns.py`
- `tool_comparison/scripts/add_sequence_accuracy_columns.py` (not currently invoked)
- `tool_comparison/hail_batch_pipelines/atarva_pipeline.py` (reference only, not invoked directly --
  `06_score_strts_native.py` reproduces its `create_atarva_step` command line)

They implement str-truth-set's own allele-length concordance comparison (the same code behind
the [tool comparison viewer](https://broadinstitute.github.io/str-truth-set-v2/tool_comparison_viewer.html)'s
rankings), reused here via `../06_score_strts_native.py` so our benchmark's "length-only" metric
matches that methodology exactly rather than approximating it.

Requires `str_analysis` (pinned commit, for the ATaRVa-VCF-to-ExpansionHunter-JSON converter
and the JSON-to-TSV combiner) and `pandas<3` (these scripts predate pandas 3's strict string
dtype and raise a `TypeError` under it) -- see `06_score_strts_native.py`'s docstring.
