# Reproducing All Tables and Figures in the Paper

This subfolder provides a single entry point to reproduce every table and figure in:

> **GenoAdjust: genome-scale protein function adjustment using constraint optimization**  
> Asaad Mohammedsaleh and Robert Hoehndorf (*Bioinformatics*, `BIOINF-2026-2454`)

---

## 1. Instant Regeneration of Tables from Committed Summary Artifacts (< 5 seconds, CPU-only)

All aggregated JSON/TSV evaluation outputs from our 564-proteome experiments are committed under `analysis/revision_results/` and `analysis/revision_results/rev2/`. To regenerate formatted Markdown, CSV/TSV, and LaTeX tables into `reproducibility/output/`:

```bash
uv run python reproducibility/reproduce_all_tables_and_figures.py
```

---

## 2. Full End-to-End Cluster Reproduction (SLURM)

To re-run all upstream predictors, Stage 1 and Stage 2 CP-SAT optimizations, GAEF and CAFA evaluations, and all revision experiments across the 564 UniProtKB/Swiss-Prot reference proteomes on a SLURM cluster:

```bash
cp config.env.example config.env
# Edit config.env to set DATA_DIR, GAEF_DIR, DEEPGO_SE_REPO, SPROF_GO_REPO
./data_prep/download_data.sh
bash reproducibility/run_all_slurm.sh
```

---

## 3. Complete Map of Paper Tables and Figures to Code and Data

### Main Manuscript (`paper-draft.pdf`)

| Item | Description | Producing Script / SLURM Job | Committed Summary Data |
| :--- | :--- | :--- | :--- |
| **Table 1** (`tab:flips`) | Stage 1 & Stage 2 annotation flips and rates across all predictors and splits | `analysis/rev2_tau_and_flips.py` (`slurm/rev2_table1_flips.slurm`) | `analysis/revision_results/rev2/flips_*.json` $\to$ `reproducibility/output/Table3_flips.md` & `Main_Table1.tex` |
| **Table 2** (`tab:gaef_preds`) | GAEF metrics before and after Stage 1 & Stage 2 (`Standard` and `--demotion-only`) | `analysis/rev2_stage2_and_baselines.py` (`slurm/rev2_stage2_and_baselines.slurm`) | `analysis/revision_results/rev2/stage2_and_baselines_*.json` $\to$ `reproducibility/output/Main_Table2.tex` |

### Supplementary Material (`paper-supplementary.pdf`)

| Item | Description | Producing Script / SLURM Job | Committed Summary Data |
| :--- | :--- | :--- | :--- |
| **Table S1** | Full-proteome CAFA metrics on the 169-protein timeset benchmark | `pipeline/evaluate_directory.py` | `reproducibility/output/Supp_Table_S1.tex` |
| **Table S2** | Full-proteome CAFA metrics on the 219-organism heldout split | `pipeline/evaluate_directory.py` | `reproducibility/output/Supp_Table_S2.tex` |
| **Table S3** | Stage 1 CP-SAT solver runtimes across CC, MF, BP | `taxon_consistency/adjust_ortools.py` | `reproducibility/output/Supp_Table_S3.tex` |
| **Table S4** | Stage 2 CP-SAT solver runtimes on CC | `complex_coherence/adjust_ortools.py` | `reproducibility/output/Supp_Table_S4.tex` |
| **Figures S1–S2** | Metagenomic taxonomic recovery across predictors and sub-ontologies | `analysis/metagenomics_plot.py` | `paper/oup-authoring-template/figures/metagenomics_fig_*.pdf` |
| **Figure S3** | TikZ diagram of the $\mathcal{V}_{\mathrm{H2}}$ sparse cone reduction | Compiled directly in `paper-supplementary.tex` | `paper-supplementary.tex` (Section S5) |
| **Table S5** | Heuristic ablation (H1/H2/H3) on taxon `759272` | `complex_coherence/heuristic_ablation.py` (`slurm/ablation/heuristic.slurm`) | `reproducibility/output/Supp_Table_S5.tex` |
| **Table S6** | Heuristic ablation across small, medium, and large proteomes | `analysis/ablation_var_counts.py` | `analysis/revision_results/ablation_var_counts.csv` |
| **Table S7** | Top-20 most-violated taxon constraints in curated Swiss-Prot ground truth | `analysis/top_violated_pairs.py` | `analysis/revision_results/top20_violated_pairs.tsv` |
| **Table S8** | Ground-truth GAEF breakdown across 514 reference proteomes | `slurm/eval/gaef_annotations.slurm` | `reproducibility/output/Supp_Table_S8.tex` |
| **Table S9** | Per-protein CAFA metrics on modified proteins (heldout split) | `pipeline/evaluate_directory.py` | `reproducibility/output/Supp_Table_S9.tex` |
| **Table S10** | 8-point decision threshold ($\tau \in [0.001, 0.5]$) sensitivity sweep | `analysis/rev2_tau_and_flips.py` (`slurm/rev2_tau_sweep.slurm`) | `analysis/revision_results/rev2/tau_sweep_*.json` $\to$ `reproducibility/output/TableS10_tau_sweep.md` |
| **Table S11** | Flip-rate reconciliation (50 timeset vs. 564 proteomes) & cost ablation (`margin`, `uniform`, `ic`) | `analysis/flip_rate.py` + `analysis/compute_go_ic.py` | `analysis/revision_results/cost_ablation.txt` + `analysis/revision_results/tau_sweep_cc.tsv` |
| **Table S12** | 7-rank NCBI taxonomic lineage inference & cross-aspect consistency | `analysis/evaluate_lineage_inference.py` (`slurm/rev2_lineage_and_funtaxis.slurm`) | `analysis/revision_results/rev2/lineage_inference.{json,tsv}` $\to$ `reproducibility/output/TableS12_lineage_inference.md` |
| **Table S13** | Comparison vs. Naive filter, `FunTaxIS-lite`, Soft SBR ($\lambda \in \{0.2, 0.5\}$), and Stage 2 (`Standard` vs. `--demotion-only`) | `analysis/rev2_stage2_and_baselines.py` + `analysis/funtaxis_lite_filter.py` | `analysis/revision_results/rev2/stage2_and_baselines_*.json` $\to$ `reproducibility/output/TableS13_baselines_soft_stage2.md` |
| **Table S14** | Soft-regularizer ($\lambda$) sweep across 564 DeepGO-SE CC proteomes | `analysis/soft_regularizer.py` | `analysis/revision_results/soft_regularizer_deepgo-se_cc.tsv` |
| **Table S15** | Complex Portal subunit validation of Standard vs. Demotion-Only Stage 2 | `analysis/rev2_stage2_and_baselines.py` + `analysis/complex_validation/validate_stage2.py` | `analysis/revision_results/rev2/cp_val_*.json` $\to$ `reproducibility/output/TableS15_complex_portal.md` |
