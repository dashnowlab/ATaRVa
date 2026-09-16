# ATaRVa HG002 benchmark

Scores ATaRVa's genotype calls against the [str-truth-set-v2](https://github.com/broadinstitute/str-truth-set-v2)
HG002 truth set, on a small subset of loci.

## The subset

`scripts/01_select_loci.py` selects the locus catalog (`data/loci.bed`, seed 42), restricted to
loci fully contained in HG002's dipcall-confident regions. The loci are randomly sampled as follows:
- 1000 loci with 2-6 bp motifs (classic STRs)
- 100 homopolymers (1 bp motifs)
- 100 VNTRs (motif_len>=7)
1200 total loci

Each locus is padded 1500bp and overlapping windows merged into disjoint regions, each
becoming its own small "pseudo-contig" (`<chrom>s<start>e<end>` -- letters, not underscores,
as separators: ATaRVa v0.7.1's own region-file validation requires the CHROM column to be
`str.isalnum()`) in a synthetic reference. Reads (soft-clipped to the window, full original
SEQ/QUAL kept), the truth VCF, and the TR catalog are all remapped into the same coordinate
space (`scripts/02_build_subset.py`), giving a small, self-contained, portable bundle. Data sources
(str-truth-set-v2's own, so results are consistent with their conventions): hg38 reference,
`HG002.bed.gz` (catalog), pacbio HiFi at 30x / ONT at native ~26x coverage,
`HG002.tandem_repeat_genotypes.tsv.gz` / the corresponding `HG002.tandem_repeats.vcf.gz`
(TR-only truth, *not* the genome-wide small-variant VCF -- see `02_build_subset.py`'s
`TRUTH_VCF_URL` comment for why that distinction matters), and `HG002.dip.bed.gz` (confident
regions).

## Running it

```bash
# one-time setup, per ATaRVa version/venv you want to benchmark, e.g.:
python3 -m venv benchmark/atarva_env && source benchmark/atarva_env/bin/activate
pip install build && pip install .   # from the repo root; needs autoconf/automake/libtool
                                      # (brew install autoconf automake libtool) to build parasail
# -- or, for a plain PyPI release: pip install ATaRVa==<version>

# once per machine, in whichever venv you'll run scoring from:
pip install 'git+https://github.com/broadinstitute/str-analysis@3d5e3dc37161d41a1bc6f92afdcf3fc99e81b30a' \
    intervaltree 'pandas<3' pysam

# build the subset (once; ~15-20 min, mostly network-bound remote region queries)
python3 benchmark/scripts/01_select_loci.py
python3 benchmark/scripts/02_build_subset.py --tech pacbio
python3 benchmark/scripts/02_build_subset.py --tech ont --skip-shared

# run ATaRVa + score (per technology, per ATaRVa venv/version you want to compare)
benchmark/scripts/03_run_atarva.sh pacbio benchmark/data benchmark/atarva_env
# ^ prints the version it detected and the exact score commands to run next, e.g.:
python3 benchmark/scripts/04_score.py --tech pacbio --atarva-version 0.7.1+ext0.01
python3 benchmark/scripts/06_score_strts_native.py --tech pacbio --atarva-version 0.7.1+ext0.01
```

Aardvark (`04_score.py`) needs Docker -- it runs the linux-x86_64 binary under `--platform
linux/amd64` (a no-op on a native x86_64 CI runner, emulated on Apple Silicon), staging inputs
through a plain system temp dir before the bind-mount, since mounting straight from a
cloud-synced checkout (OneDrive, Dropbox, iCloud Drive) can intermittently deadlock on macOS.
`06_score_strts_native.py` needs `str_analysis` and `pandas<3` in the same venv as ATaRVa (see
above).

## Comparing ATaRVa versions over time

`03_run_atarva.sh` writes a **version-tagged** output (`data/atarva.<tech>.<version>.vcf.gz`,
version auto-detected from `atarva[-ext] --version`), so multiple ATaRVa builds' outputs can
coexist in `data/` and be rescored any time without re-running genotyping. Every scoring run
(`04_score.py` and `06_score_strts_native.py`, both requiring `--atarva-version`) *appends* its
metrics as rows to `baseline_history.tsv` (`scripts/baseline_history.py`) rather than
overwriting a single snapshot file -- so results from every version ever benchmarked stay
around. Filter/pivot it to compare versions, e.g.:

```bash
awk -F'\t' '$3=="pacbio" && $6=="BASEPAIR/f1"' benchmark/baseline_history.tsv
```

**CI only ever benchmarks the version it's checking out** (a release tag, or `main` on the
monthly schedule) -- see "Continuous integration" below. Comparing against a *different* ATaRVa
version (an upstream release, an older tag) is a local-only activity: build/install that other
version into its own venv (a fresh `atarva_env_*` directory works well -- see
`atarva_env_v0.7.1/` for the pattern, gitignored like `atarva_env/`), run `03_run_atarva.sh`
with that venv for both technologies, then `04_score.py` and `06_score_strts_native.py` with the
version it printed -- then commit the resulting `baseline_history.tsv` rows in a normal PR if
you want that comparison kept. The subset itself (loci, reference, reads, truth) doesn't need
rebuilding between versions -- only re-genotyping and re-scoring.

## Continuous integration

`.github/workflows/benchmark.yml` runs this benchmark on `ubuntu-latest`, triggered by a
published release, a monthly schedule (the 1st, 06:00 UTC), or manually
(`workflow_dispatch`). Per run:

1. Restores `data/` and `research/` from `actions/cache` (keyed on a hash of
  `01_select_loci.py` + `02_build_subset.py`, so a deliberate catalog-design change busts the
  cache automatically); builds them fresh only on a cache miss (~15-20 min, otherwise skipped).
2. Builds ATaRVa from the checked-out source into a fresh venv, genotypes both technologies,
  scores both ways.
3. Runs `07_check_regression.py` against the checked-in `baseline_current.tsv` -- **fails the
  job** if any guarded metric dropped by more than its tolerance (0.01 by default; see that
  script's docstring for why the tolerance is looser than str-truth-set-v2's own 0.002 -- our
  n=1196 is much smaller than their n=146k+, so a rerun has more room for legitimate noise).
4. Uploads `baseline_history.tsv` and `data/latest_metrics.tsv` as a build artifact.

It does **not** auto-commit results back to the repo (would need CI write access and could race
concurrent runs). To keep a CI run's numbers in history, download the artifact and commit
`baseline_history.tsv` by hand. To accept a genuine, intentional change in the accepted
baseline, run `07_check_regression.py --update-baseline` locally and commit
`baseline_current.tsv` with a message explaining why -- same discipline as
str-truth-set-v2's own `truthset_release_check_baseline.tsv`.

## Two scoring methodologies

- **`aardvark_basepair`** (`04_score.py`): full-**sequence** edit-distance comparison
  (Aardvark's BASEPAIR metric only -- its stricter, hap.py-style GT metric is deliberately
  excluded; see `guarded_metrics()`'s comment for why it isn't a useful signal for TR loci).
  This reuses the approach of str-truth-set-v2's `docs/truthset_release_check.md` -- but that
  doc uses it only to QA the truth set itself against an independent GIAB benchmark, not to
  rank tools.
- **`strts_native`** (`06_score_strts_native.py`): **allele-length** (repeat copy number)
  concordance, reusing str-truth-set's own comparison scripts (`scripts/strts_reused/`, see its
  README) unmodified. This is what actually produces the [tool comparison viewer](https://broadinstitute.github.io/str-truth-set-v2/tool_comparison_viewer.html)'s
  published ATaRVa numbers. Reports both a strict exact-match rate and a +/-1-repeat-unit
  tolerant rate (the ATaRVa paper's own headline metrics are similarly exact-match and
  +/-1bp-tolerant).

Sequence match is a strictly harder bar than length match (a correct repeat count with a
different internal motif arrangement still counts as correct length-wise but not
sequence-wise), so `aardvark_basepair` numbers are expected to run lower than `strts_native`'s
for the same run. Neither comparison is polluted by nearby non-repetitive variants (e.g. local
SNPs): both truth (`HG002.tandem_repeats.vcf.gz`) and ATaRVa's own output are TR-only VCFs, so
there's nothing else for either comparison to pick up regardless of how much flanking sequence
is included in a scored region.

Both scripts default `--catalog` to `random_2to6bp_with_extremes`; pass a different label only
if you build a different catalog design again in the future.

## Current results (2026-09-14)

`random_2to6bp_with_extremes`, n=1196 (scored) of 1200 (selected) -- all motif lengths and
zygosities mixed (not restricted to heterozygous truth loci, unlike some published headline
figures). ATaRVa is run with `--karyotype XY --min-reads 2` (not its own defaults), matching
str-truth-set-v2's `atarva_pipeline.py` exactly -- HG002 is male, so chrX outside the PAR and
chrY need haploid genotyping, and the lower `--min-reads` avoids unfairly no-calling
marginal-coverage loci relative to how other tools are benchmarked there.

| metric | ONT 0.7.1 | ONT 0.7.1+ext0.01 | PacBio 0.7.1 | PacBio 0.7.1+ext0.01 |
|---|---|---|---|---|
| BASEPAIR f1 (Aardvark, full sequence) | 0.85 | 0.78 | 0.77 | 0.79 |
| BASEPAIR precision | 0.92 | 0.92 | 0.94 | 0.94 |
| BASEPAIR recall | 0.79 | 0.68 | 0.65 | 0.68 |
| ExactMatch, all 1196 loci | 42.6% | 44.5% | 82.0% | 84.1% |
| ExactMatch, 2-6bp only (n=1000) | 45.3% | 46.9% | 83.0% | 85.4% |
| LengthMatch +/-1 unit, all 1196 | 80.7% | 80.9% | 94.8% | 95.2% |
| LengthMatch +/-1 unit, 2-6bp only | 85.2% | 85.4% | 94.7% | 95.1% |

The fork (0.7.1+ext0.01) and upstream (0.7.1) track each other closely on PacBio. On ONT,
upstream's BASEPAIR F1 (0.85) is notably ahead of the fork's (0.78) despite the fork's
length-based metrics being marginally *better* -- worth a closer look at specific loci
(`benchmark/data/strts_comparison/reads.ont.<version>.with_concordance.tsv`) before treating it
as a confirmed regression, since BASEPAIR and length-match disagreeing in direction usually
means a representation difference rather than a wrong genotype.

`ExactMatch/ALL` sits consistently below `ExactMatch/2to6bp` here (by ~1-3 points) -- the 200
homopolymer+VNTR loci (17% of the catalog) are measurably harder than the 2-6bp core, as
expected, without swamping the aggregate the way an unstratified genome-wide sample would (that
was the point of adding them deliberately rather than sampling naturally).

Cross-checked against two independent published sources: the tool comparison viewer's ATaRVa
ONT/26x/2-6bp/HET-only exact-match figure is 70.6% (146,672 loci) -- our 2-6bp-only figure of
~46-47% and the ATaRVa paper's own PacBio HiFi 30x exact-match (86.4%, close to our 83-85%) both
land in a plausible range once compared like-for-like (motif bucket), with the residual gap
explained by sample size (1000 vs 146k+ loci) and zygosity mix (not HET-restricted here), not a
pipeline bug.

## Known limitations

- `data/` (the built bundle: reference, BAMs, VCFs) is gitignored and regenerated locally or by
  CI's `actions/cache` -- ~900MB, governed mostly by the two BAMs, well past a git-friendly
  size, because reads keep their full original SEQ/QUAL rather than being trimmed to the
  window.
- Not restricted to heterozygous truth loci (unlike some published headline figures) --
  deliberately, since a real benchmark shouldn't cherry-pick the easier zygosity class.
- The catalog is a single draw (seed 42) rather than a repeated/bootstrapped sample, so treat
  percentage-point differences as indicative, not statistically confirmed.
- CI benchmarks only the version it checks out -- no automated regression testing against other
  ATaRVa versions (upstream releases, older tags); those comparisons are local-only (see
  "Comparing ATaRVa versions over time").
