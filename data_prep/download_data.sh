#!/usr/bin/env bash
# Fetch large constraint files that are not committed to the repo.
#
# Usage:
#   cp config.env.example config.env       # then edit
#   source config.env
#   ./data_prep/download_data.sh
#
# Pinned releases (update intentionally; do not silently track latest):
#   go-basic.obo: 2025-10 release (matches paper)
#   ncbitaxon.obo: 2025-09 release (matches paper)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[ -f "${REPO_ROOT}/config.env" ] && source "${REPO_ROOT}/config.env"

DATA_LOCAL="${REPO_ROOT}/data"
mkdir -p "${DATA_LOCAL}"

GO_OBO_URL="${GO_OBO_URL:-http://release.geneontology.org/2025-10-01/ontology/go-basic.obo}"
NCBITAXON_OBO_URL="${NCBITAXON_OBO_URL:-http://purl.obolibrary.org/obo/ncbitaxon.obo}"

if [ ! -f "${DATA_LOCAL}/go-basic.obo" ]; then
    echo "==> Fetching go-basic.obo"
    curl -L --fail -o "${DATA_LOCAL}/go-basic.obo" "${GO_OBO_URL}"
else
    echo "==> go-basic.obo already present, skipping"
fi

if [ ! -f "${DATA_LOCAL}/ncbitaxon.obo" ]; then
    echo "==> Fetching ncbitaxon.obo (large, ~700 MB)"
    curl -L --fail -o "${DATA_LOCAL}/ncbitaxon.obo" "${NCBITAXON_OBO_URL}"
else
    echo "==> ncbitaxon.obo already present, skipping"
fi

if [ ! -f "${DATA_LOCAL}/ncbitaxon_hierarchy.tsv" ]; then
    echo "==> Parsing ncbitaxon.obo -> ncbitaxon_hierarchy.tsv"
    python "${REPO_ROOT}/data_prep/parse_ncbitaxon.py" \
        --ncbitaxon-obo-file "${DATA_LOCAL}/ncbitaxon.obo" \
        --output-file "${DATA_LOCAL}/ncbitaxon_hierarchy.tsv"
else
    echo "==> ncbitaxon_hierarchy.tsv already present, skipping"
fi

echo ""
echo "Done. Bundled + downloaded data in ${DATA_LOCAL}/:"
ls -lh "${DATA_LOCAL}/"

cat <<'EOF'

Next steps:
  - The UniProt reference proteomes and Swiss-Prot annotations are not fetched
    by this script. The repo ships `data/uniprot_proteomes_ids.tsv`; rebuild it
    only if you re-derive the benchmark. To rebuild the per-taxon annotations
    from already-downloaded reference proteomes:
      python data_prep/set_up_proteomes_dir.py \
          --main_proteomes_dir   ${DATA_DIR}/uniprot_reference_proteomes \
          --output_dir           ${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic \
          --go_file              data/go-basic.obo \
          --proteomes_ids_file   data/uniprot_proteomes_ids.tsv
    Expect tens of GB and several hours of disk I/O.
  - For the baselines that need pretrained weights:
      DeepGO-SE:  see baselines/deepgo_se/README.md
                  (run inference via slurm/train/predict_deepgo_se.slurm)
      SPROF-GO:   see baselines/sprof_go/README.md
                  (run inference via slurm/train/predict_sprof_go.slurm)
      MLP-ESM2:   train + predict via slurm/train/train_mlp_heldout.slurm
                  (ESM2 embeddings produced first by slurm/train/extract_esm2.slurm)
EOF
