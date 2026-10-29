# genome-scale-pfp-adjust

Post-hoc CP-SAT constraint optimization for repairing genome-scale protein function predictions.

Companion code for the paper *Genome-scale protein function adjustment using constraint optimization*. Implements a two-stage CP-SAT solver (taxon consistency, then complex coherence) that seeks low-cost edits of genome-scale GO-term predictions under explicitly encoded constraints.

## What the solver does

Given per-protein GO predictions from any upstream method, the solver seeks a low-cost set of edits satisfying its encoded constraints:

1. **Stage 1 — taxon consistency.** No GO term in the predicted set violates the organism's NCBI Taxonomy lineage (e.g. no *photosynthesis* in an animal proteome). Inferred from `only_in_taxon` / `never_in_taxon` constraints in the Gene Ontology.

2. **Stage 2 — complex coherence.** No term in the operational complex vocabulary is annotated to exactly one protein in the proteome. Singletons are either demoted from the lone protein or promoted onto another protein. This does not establish that the proteins are interaction partners.

Both stages are encoded as integer programs and solved with OR-Tools CP-SAT. Stage 2 closes a sparse repair domain over hierarchy and affected complex constraints. Top-k candidate restriction is a heuristic. Bounded runs distinguish feasible solutions from proven optima and record objective bounds.

An optional conservative mode restricts edits to supplied evidence, preserves
protected annotations, and reports unresolved singletons. Unknown-taxon solves
can also report which ranked taxa are required across every optimal GO repair.
The [extension protocol](reproducibility/EXTENSIONS.md) describes these options,
the standalone repair command, and the controlled recovery evaluation.
The [validated extension results](analysis/revision_results/extensions/README.md)
report preservation, all-optimum lineage conclusions and the recovery controls.

## Demo

An executable Showboat demo is available in [`demos/empty_quarter_ibex/README.md`](demos/empty_quarter_ibex/README.md). It fetches a real Empty Quarter genome annotation from the public reproducibility bundle at `https://bio2vec.net/data/genoadjust/rh04_63_rh04_bacillus_spizizenii_genoadjust_demo/` (`rh04` / `63_rh04`, *Bacillus spizizenii*, `NCBITaxon_96241`), converts PGAP GO annotations into the solver's prediction format, and verifies both taxon-consistency and complex-coherence repairs.

## External dependencies

Two pretrained baselines and the evaluator use code from external repositories. Clone them and point `config.env` at the clones:

- `GAEF_DIR` — `https://github.com/bio-ontology-research-group/GAEF` (the current evaluator imports `GAEF.completeness` from upstream; the divergent modules — `coherence`, `taxon_consistency`, `complex_classifier`, and the `Ontology` IC methods — live under `gaef_patches/` in this repo, and the constraint files under `data/constraints/`).
- `DEEPGO_SE_REPO` — `https://github.com/bio-ontology-research-group/deepgo2`
- `SPROF_GO_REPO` — `https://github.com/biomed-AI/SPROF-GO`

## Quickstart

```
# 1. CPU revision environment (upstream retraining is separate)
python3.10 -m venv /your/project/revision-venv
source /your/project/revision-venv/bin/activate
pip install -r requirements-revision.txt

# 2. Local config
cp config.env.example config.env
$EDITOR config.env                              # set DATA_DIR, CONDA_*, GAEF_DIR, …
source config.env

# 3. Fetch large constraint files
./data_prep/download_data.sh                    # go-basic.obo, ncbitaxon.obo

# 4. Tests (offline, CPU-only, <30 s) — also run in GitHub Actions CI
python data/constraints/verify_checksums.py   # pinned constraint files
python tests/test_data_loaders.py
python tests/test_solver_smoke.py
python tests/test_golden_e2e.py               # 50-protein organism, diffs vs golden output
# (these four need only `pip install -r requirements-test.txt`, not the full conda env)

# 5. Adjustment/evaluation on fixed predictor inputs
# Configure reproducibility/revision.env.example, then:
bash reproducibility/run_all_slurm.sh /absolute/path/to/revision.env

# Regenerate local Seq-Sim and MLP inputs first using reproducibility/BASELINES.md.
# Historical local-baseline training wrappers are retired.

```

## Repository layout

```
taxon_consistency/        Stage 1 CP-SAT solver
complex_coherence/        Stage 2 CP-SAT solver + heuristic ablation
pipeline/                 Adjustment orchestrator (chains both stages) + evaluation scripts
gaef_patches/             Local extensions to upstream GAEF (see "External dependencies")
deepgo/                   Utility package (Ontology, MLPBlock, …)
baselines/
  seq_sim/                Sequence-similarity baseline (DIAMOND BLAST)
  mlp_esm2/               MLP trained on ESM2 embeddings.
  deepgo_se/              Wrapper for DeepGO-SE (external repo)
  sprof_go/               Wrapper for SPROF-GO (external repo)
data_prep/                Dataset/split construction, OBO parsers, dataset builder, download_data.sh
data/                     Bundled small constraint TSVs (large ones downloaded)
  constraints/            GAEF constraint files (vendored — see gaef_patches/)
splits/{timeset,heldout}/ Train/test organism splits used in the paper
analysis/                 Supplementary & revision analyses (term-centric, metagenomics, flip ranking, baselines)
reproducibility/          Single-command reproduction of all paper & supplementary tables and figures
slurm/{adjust,eval,train,ablation}/
                          Cluster job scripts, all templated via config.env
tests/                    CPU-only smoke tests + fixtures
```

## `DATA_DIR` layout

`DATA_DIR` is the single working-data root referenced by every SLURM script and
CLI default. You populate `uniprot_reference_proteomes/` yourself (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/); everything under `swissprot_proteomes_folds/` is produced by
`data_prep/`, the training/prediction jobs, and the adjustment stages.

```
${DATA_DIR}/
├── uniprot_reference_proteomes/        # UniProt reference proteomes as downloaded from the UniProt FTP site
└── swissprot_proteomes_folds/          # all generated artifacts live here
    ├── annotations-go-basic/           # ground-truth GO annotations (data_prep/)
    ├── esm2_embeddings/                # per-protein ESM2 mean embeddings (train/extract_esm2.slurm)
    ├── seq_sim_{timeset,heldout}_results/{cc,mf,bp}/
    │   ├── predictions/                # raw upstream predictions
    │   ├── optimized/                  # after Stage 1 (taxon) + Stage 2 (complex, cc only)
    │   └── optimized_without_taxons/   # after Stage 1 and 2, but without providing taxon information
    ├── mlp_heldout_results/{cc,mf,bp}/<same subdirs>
    ├── deepgo-se_results/{cc,mf,bp}/<same subdirs>
    └── sprof-go_results/{cc,mf,bp}/<same subdirs>
```

SPROF-GO additionally uses `${DATA_DIR}/sprof_go/sprof_go_work/` as an
intermediate working directory (configurable via `WORK_DIR=` on the sbatch
line).

## Reproducing the revised analysis

The October 2026 revision repairs score export, taxonomy-union semantics, evaluation membership and accuracy metrics. Previous numerical tables require recomputation. See [the revision workflow](reproducibility/README.md) for required inputs, the complete task manifest and the status of the result bundle.

```bash
# Adjustment and evaluation on fixed predictor inputs; configure paths first.
bash reproducibility/run_all_slurm.sh /absolute/path/to/revision.env

# Generate tables only after the complete numerical bundle is available.
python reproducibility/reproduce_all_tables_and_figures.py \
  --results /path/to/complete-results --output /path/to/generated
```

The generator rejects missing results and does not extract numbers from the manuscript. Regenerating the local Seq-Sim and MLP baselines is a prerequisite described in [BASELINES.md](reproducibility/BASELINES.md). Supplied DeepGO-SE and SPROF-GO predictions remain fixed. Public archival verification is a separate final step.

## Citation

> Mohammedsaleh, A. and Hoehndorf, R. (2026). Genome-scale protein function adjustment using constraint optimization. (under preparation).
