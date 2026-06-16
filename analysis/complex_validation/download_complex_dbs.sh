#!/usr/bin/env bash
# Fetch curated complex subunit lists used to validate Stage 2 promotions
# and demotions. Run from this directory or pass DATA_DIR via config.env.
#
# Usage:
#   source ../../config.env
#   bash download_complex_dbs.sh
#
# Outputs into ${DATA_DIR}/complex_validation/:
#   complexportal/<taxon>.tsv      -- ComplexPortal complextab TSVs
#   corum/coreComplexes.txt        -- CORUM core complexes
#
# Note: CORUM may require updating the URL if the maintainers move the file.
# If the curl below fails, fetch coreComplexes.txt manually from
#   https://mips.helmholtz-muenchen.de/corum/
# and drop it into corum/.

set -euo pipefail

if [ -z "${DATA_DIR:-}" ]; then
    echo "ERROR: DATA_DIR is not set. Source config.env first." >&2
    exit 1
fi

DEST="${DATA_DIR}/complex_validation"
mkdir -p "${DEST}/complexportal" "${DEST}/corum"

# ---------------------------------------------------------------------------
# EBI Complex Portal — one TSV per taxon
# ---------------------------------------------------------------------------
# Index: https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/
# The taxon list below covers the model + reference organisms most likely to
# overlap with the 50 timeset organisms (and additional well-curated species).
# Append more taxon IDs here if you want broader coverage.
CP_BASE="https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab"
CP_TAXA=(
    83333    # Escherichia coli K-12
    559292   # Saccharomyces cerevisiae S288C
    284812   # Schizosaccharomyces pombe 972h-
    6239     # Caenorhabditis elegans
    7227     # Drosophila melanogaster
    7955     # Danio rerio
    8355     # Xenopus laevis
    9606     # Homo sapiens
    10090    # Mus musculus
    10116    # Rattus norvegicus
    3702     # Arabidopsis thaliana
    4577     # Zea mays
    224308   # Bacillus subtilis 168
    1773     # Mycobacterium tuberculosis
    272634   # Mycoplasma pneumoniae
)

echo "==> Fetching Complex Portal data"
for taxon in "${CP_TAXA[@]}"; do
    out="${DEST}/complexportal/${taxon}.tsv"
    if [ -f "${out}" ]; then
        echo "    [skip] ${taxon}.tsv already present"
        continue
    fi
    url="${CP_BASE}/${taxon}.tsv"
    if curl -sSfL "${url}" -o "${out}"; then
        echo "    [ok]   ${taxon}.tsv"
    else
        echo "    [miss] no ComplexPortal data for taxon ${taxon}"
        rm -f "${out}"
    fi
done

# ---------------------------------------------------------------------------
# CORUM core complexes
# ---------------------------------------------------------------------------
echo "==> Fetching CORUM coreComplexes"
CORUM_URL="${CORUM_URL:-https://mips.helmholtz-muenchen.de/corum/download/releases/current/coreComplexes.txt.zip}"
if [ ! -f "${DEST}/corum/coreComplexes.txt" ]; then
    tmp="${DEST}/corum/coreComplexes.txt.zip"
    if curl -sSfL "${CORUM_URL}" -o "${tmp}"; then
        unzip -o -q "${tmp}" -d "${DEST}/corum/"
        rm -f "${tmp}"
        echo "    [ok]   coreComplexes.txt"
    else
        cat >&2 <<EOF
    [miss] CORUM download failed (${CORUM_URL}).
           CORUM occasionally moves the download URL or requires registration;
           fetch coreComplexes.txt manually from https://mips.helmholtz-muenchen.de/corum/
           and place it at ${DEST}/corum/coreComplexes.txt.
EOF
    fi
else
    echo "    [skip] coreComplexes.txt already present"
fi

echo
echo "Done. Inventory in ${DEST}/:"
ls -lh "${DEST}/complexportal/" "${DEST}/corum/" 2>/dev/null || true
