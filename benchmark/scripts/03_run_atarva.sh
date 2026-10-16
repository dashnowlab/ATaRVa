#!/usr/bin/env bash
# Run ATaRVa's genotype command against the benchmark subset for one technology, writing a
# version-tagged output file so multiple ATaRVa versions' results can coexist and be
# rescored/compared without re-running genotyping (see baseline_history.py).
#
# Works against any ATaRVa install with the -ext fork's `atarva-ext` command or upstream's
# plain `atarva` command -- whichever is on PATH after activating atarva-env.
#
# Usage: 03_run_atarva.sh <pacbio|ont> [data-dir] [atarva-env]
set -euo pipefail

TECH="${1:?usage: 03_run_atarva.sh <pacbio|ont> [data-dir] [atarva-env]}"
DATA_DIR="${2:-benchmark/data}"
ATARVA_ENV="${3:-benchmark/atarva_env}"

source "${ATARVA_ENV}/bin/activate"

if command -v atarva-ext >/dev/null 2>&1; then
  ATARVA_BIN=atarva-ext
elif command -v atarva >/dev/null 2>&1; then
  ATARVA_BIN=atarva
else
  echo "No atarva-ext or atarva command found on PATH after activating ${ATARVA_ENV}" >&2
  exit 1
fi

VERSION="$("${ATARVA_BIN}" --version 2>&1 | tail -1 | sed -E 's/^ATaRVa version //')"
echo "Using ${ATARVA_BIN} ${VERSION}"

OUT_PREFIX="${DATA_DIR}/atarva.${TECH}.${VERSION}"

# --karyotype XY: HG002 is male, so chrX outside the PAR and all of chrY must be genotyped
# haploid, matching str-truth-set-v2's own atarva_pipeline.py -- without this, ATaRVa treats
# every locus as diploid and every hemizygous truth call becomes an unmatchable false call.
# --min-reads 2 (not ATaRVa's own default of 10): also matches atarva_pipeline.py, which
# lowers this so low-coverage/haploid loci aren't unfairly no-called relative to other tools
# benchmarked at a lower threshold.
"${ATARVA_BIN}" genotype \
  -f "${DATA_DIR}/ref.fa" \
  -b "${DATA_DIR}/reads.${TECH}.bam" \
  -r "${DATA_DIR}/catalog.bed.gz" \
  -o "${OUT_PREFIX}.vcf" \
  --karyotype XY \
  --min-reads 2 \
  -t 4

# ATaRVa's own output isn't always coordinate-sorted across our many small pseudo-contigs
# (seen concretely with v0.7.1 on the ONT run) -- sort defensively before indexing rather
# than assume it.
{ grep "^#" "${OUT_PREFIX}.vcf"; grep -v "^#" "${OUT_PREFIX}.vcf" | sort -k1,1 -k2,2n; } \
  > "${OUT_PREFIX}.sorted.vcf"
mv "${OUT_PREFIX}.sorted.vcf" "${OUT_PREFIX}.vcf"

bgzip -f "${OUT_PREFIX}.vcf"
tabix -p vcf "${OUT_PREFIX}.vcf.gz"

echo "Wrote ${OUT_PREFIX}.vcf.gz"
echo "Score it with: python3 benchmark/scripts/04_score.py --tech ${TECH} --atarva-version ${VERSION}"
echo "           and: python3 benchmark/scripts/06_score_strts_native.py --tech ${TECH} --atarva-version ${VERSION}"
