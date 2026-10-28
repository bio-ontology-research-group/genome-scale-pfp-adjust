#!/usr/bin/env bash
# =============================================================================
# Master SLURM Orchestration Script to Reproduce All Experiments, Tables, and
# Figures in GenoAdjust (Bioinformatics, BIOINF-2026-2454)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

if [[ ! -f config.env ]]; then
    echo "ERROR: config.env not found. Copy config.env.example to config.env and set DATA_DIR." >&2
    exit 1
fi
source config.env

echo "=== Step 1: Verify pinned constraint checksums ==="
python data/constraints/verify_checksums.py

echo "=== Step 2: Upstream Predictions (ESM2, Seq-Sim, MLP-ESM2, DeepGO-SE, SPROF-GO) ==="
ESM_JOB=$(sbatch --parsable slurm/train/extract_esm2.slurm)
SS_H_JOB=$(sbatch --parsable slurm/train/predict_seq_sim_heldout.slurm)
SS_T_JOB=$(sbatch --parsable slurm/train/predict_seq_sim_timeset.slurm)
for T in cc mf bp; do
    sbatch --dependency=afterok:${ESM_JOB} slurm/train/train_mlp_heldout.slurm "${T}"
    sbatch --dependency=afterok:${ESM_JOB} slurm/train/predict_deepgo_se.slurm "${T}"
done
SPROF_JOB=$(sbatch --dependency=afterok:${ESM_JOB} --parsable slurm/train/predict_sprof_go.slurm)

echo "=== Step 3: Two-Stage CP-SAT Adjustment (Stage 1 Taxon + Stage 2 Complex Coherence) ==="
ADJ_SST=$(sbatch --dependency=afterok:${SS_T_JOB} --parsable slurm/adjust/seq_sim_timeset.slurm)
ADJ_SSH=$(sbatch --dependency=afterok:${SS_H_JOB} --parsable slurm/adjust/seq_sim_heldout.slurm)
ADJ_MLP=$(sbatch --parsable slurm/adjust/mlp_heldout.slurm)
ADJ_DGO=$(sbatch --parsable slurm/adjust/deepgo_se.slurm)
for T in cc mf bp; do
    sbatch --dependency=afterok:${SPROF_JOB} slurm/adjust/sprof_go.slurm "${T}"
done

echo "=== Step 4: GAEF & Per-Protein CAFA Evaluations (Tables 1, 2, 4; Supp. Tables S1-S9) ==="
sbatch slurm/eval/gaef_annotations.slurm
sbatch --dependency=afterok:${ADJ_SST}:${ADJ_DGO} slurm/eval/gaef_timeset_methods.slurm
for TX in 272844 759272 10116; do
    sbatch slurm/ablation/heuristic.slurm "${TX}"
done

echo "=== Step 5: Revision Experiments (Table 3; Supp. Tables S10-S15) ==="
# Table 3: Stage 1 + Stage 2 flip counts & rates across all 15 (model x aspect) pairs
sbatch --dependency=afterok:${ADJ_SST}:${ADJ_SSH}:${ADJ_MLP}:${ADJ_DGO} slurm/rev2_table1_flips.slurm

# Supplementary Table S10: 8-point decision threshold (tau) sensitivity sweep
sbatch --dependency=afterok:${ADJ_SST}:${ADJ_DGO} slurm/rev2_tau_sweep.slurm

# Supplementary Table S12 & S13: 7-rank taxonomic lineage inference + FunTaxIS-lite
sbatch --dependency=afterok:${ADJ_SST}:${ADJ_DGO} slurm/rev2_lineage_and_funtaxis.slurm

# Supplementary Tables S13 & S15: Stage 2 (Standard vs Demotion-Only), Naive, FunTaxIS-lite, Soft SBR, Complex Portal
sbatch --dependency=afterok:${ADJ_SST}:${ADJ_DGO} slurm/rev2_stage2_and_baselines.slurm

echo "All SLURM jobs submitted. Once completed, run:"
echo "  uv run python reproducibility/reproduce_all_tables_and_figures.py"
