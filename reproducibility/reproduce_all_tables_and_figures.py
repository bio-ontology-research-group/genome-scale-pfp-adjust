#!/usr/bin/env python3
"""
Reproduce all tables (Tables 1-2, Supplementary Tables S1-S15) and figures
(Supplementary Figures S1-S2) for:

  GenoAdjust: genome-scale protein function adjustment using constraint optimization
  Asaad Mohammedsaleh and Robert Hoehndorf (Bioinformatics, BIOINF-2026-2454)

Usage:
  uv run python reproducibility/reproduce_all_tables_and_figures.py
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
REV_DIR = REPO_ROOT / "analysis" / "revision_results"
REV2_DIR = REV_DIR / "rev2"
OUT_DIR = Path(__file__).resolve().parent / "output"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_artifact(name: str, content: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / name
    out_path.write_text(content, encoding="utf-8")
    return out_path


def build_table1_flips() -> str:
    """Reproduce Main Paper Table 1 (tab:flips): Stage 1 annotation flips and rates across all models/splits."""
    models = [
        ("seq_sim_timeset", "Seq-Sim (timeset, n=50)"),
        ("deepgo-se", "DeepGO-SE* (timeset, n=50)"),
        ("sprof-go", "SPROF-GO* (timeset, n=50)"),
        ("seq_sim_heldout", "Seq-Sim (heldout, n=219)"),
        ("mlp_heldout", "MLP-ESM2 (heldout, n=219)"),
    ]
    lines = [
        "# Main Paper Table 1 (tab:flips): Stage 1 Annotation Flips and Rates",
        "",
        "| Model (Split) | Aspect | Mode | Stage 1 Flips | Preds > tau | Act.% | Grid% |",
        "| :--- | :--- | :--- | ---: | ---: | ---: | ---: |",
    ]
    for m_key, m_label in models:
        for sub in ["cc", "mf", "bp"]:
            fpath = REV2_DIR / f"flips_{m_key}_{sub}.json"
            if not fpath.exists():
                continue
            rows = load_json(fpath)
            by_mode = {r["mode"]: r for r in rows}
            for mode_key, mode_label in [("no_taxons", "no taxons"), ("with_taxons", "taxons")]:
                if mode_key not in by_mode:
                    continue
                r = by_mode[mode_key]
                lines.append(
                    f"| {m_label} | {sub.upper()} | {mode_label} | {r['flips_above_tau']:,} | "
                    f"{r['preds_above_tau']:,} | {r['active_flip_rate_pct']:.2f}% | "
                    f"{r['grid_flip_rate_pct']:.4f}% |"
                )
    return "\n".join(lines) + "\n"


def build_table_s10_tau_sweep() -> str:
    """Reproduce Supplementary Table S10 (8-point threshold sensitivity sweep)."""
    lines = [
        "# Supplementary Table S10: 8-Point Decision Threshold (tau) Sensitivity Sweep",
        "",
        "| Predictor | Aspect | tau | Preds > tau | Stage 1 Flips | Act.% | Grid% | TC | Fmax | Smin | AUC | AUPR |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for m_key, m_label in [
        ("seq_sim_timeset", "Seq-Sim"),
        ("deepgo-se", "DeepGO-SE"),
        ("sprof-go", "SPROF-GO"),
    ]:
        for sub in ["cc", "mf"]:
            p = REV2_DIR / f"tau_sweep_{m_key}_{sub}.json"
            if not p.exists():
                continue
            rows = load_json(p)
            for r in rows:
                cafa = r.get("pooled_cafa") or {}
                lines.append(
                    f"| {m_label} | {sub.upper()} | {float(r['tau']):.3f} | {r['preds_above_tau']:,} | "
                    f"{r['flips_above_tau']:,} | {r['active_flip_rate_pct']:.2f}% | "
                    f"{r['grid_flip_rate_pct']:.4f}% | {r.get('gaef_taxon_consistency_after', 1.0):.2f} | "
                    f"{cafa.get('fmax', 0.0):.4f} | {cafa.get('smin', 0.0):.2f} | "
                    f"{cafa.get('avg_auc', 0.0):.4f} | {cafa.get('aupr', 0.0):.4f} |"
                )
    return "\n".join(lines) + "\n"


def build_table_s12_lineage() -> str:
    """Reproduce Supplementary Table S12 (Multi-rank taxonomic lineage inference)."""
    p = REV2_DIR / "lineage_inference.json"
    d = load_json(p)
    lines = [
        "# Supplementary Table S12: Multi-Rank Taxonomic Lineage Inference (50 Timeset Proteomes)",
        "",
        "| Predictor | Aspect | Deepest Correct (None/Dom/Phy/Cls/Ord+) | First Divergence (None/Dom/Phy/Cls+) | Exact Prefix | Cross-Aspect Compat. (Forced / Solver) |",
        "| :--- | :--- | :--- | :--- | ---: | :--- |",
    ]
    for m_key, m_label in [
        ("seq_sim_timeset", "Seq-Sim"),
        ("deepgo-se", "DeepGO-SE"),
        ("sprof-go", "SPROF-GO"),
    ]:
        md = d[m_key]
        for aspect in ["cc", "mf", "bp", "joint_all_3"]:
            ad = md[aspect]
            dc = ad["solver_deepest_correct_counts"]
            fd = ad["solver_first_divergence_counts"]
            ord_plus = sum(dc.get(r, 0) for r in ["order", "family", "genus"])
            cls_plus = sum(fd.get(r, 0) for r in ["class", "order", "family", "genus"])
            dc_str = f"{dc.get('none',0)} / {dc.get('domain',0)} / {dc.get('phylum',0)} / {dc.get('class',0)} / {ord_plus}"
            fd_str = f"{fd.get('none',0)} / {fd.get('domain',0)} / {fd.get('phylum',0)} / {cls_plus}"
            ex = ad["solver_exact_prefix_count"]
            asp_label = "Joint (All)" if aspect == "joint_all_3" else aspect.upper()
            compat_str = (
                f"{ad.get('cross_aspect_forced_compatible_count', '--')} / {ad.get('cross_aspect_solver_compatible_count', '--')}"
                if aspect == "joint_all_3"
                else "--"
            )
            lines.append(f"| {m_label} | {asp_label} | {dc_str} | {fd_str} | {ex} / 50 | {compat_str} |")
    return "\n".join(lines) + "\n"


def build_table_s13_baselines_and_stage2() -> str:
    """Reproduce Supplementary Table S13 (Baselines, FunTaxIS-lite, Soft SBR, Stage 2 modes)."""
    lines = [
        "# Supplementary Table S13: Baselines, Soft Regularizer (SBR), and Stage 2 Modes on CC",
        "",
        "| Predictor | Method / Mode | TC (n=50) | CCoh% (n=50) | TC (n=41) | CCoh% (n=41) | Fmax | Smin | AUC | AUPR | wFmax |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    spec_order = [
        ("before", "Raw (preds)"),
        ("naive_noprop", "Naive (no prop, taxons)"),
        ("naive_prop", "Naive (prop, taxons)"),
        ("funtaxis_lite_prop", "FunTaxIS-lite (prop, taxons)"),
        ("soft_s1_lam0.2", "Soft SBR Stage 1 (lambda=0.2)"),
        ("soft_s2_lam0.2", "Soft SBR Stage 2 (lambda=0.2)"),
        ("soft_s1_lam0.5", "Soft SBR Stage 1 (lambda=0.5)"),
        ("soft_s2_lam0.5", "Soft SBR Stage 2 (lambda=0.5)"),
        ("stage1_no_taxons", "Stage 1 (no taxons)"),
        ("stage1_with_taxons", "Stage 1 (taxons)"),
        ("existing_stage2_standard_withtax", "Stage 2 Standard (taxons)"),
        ("stage2_demotion_only_withtax", "Stage 2 Demotion-Only (taxons)"),
    ]
    for m_key, m_label in [
        ("seq_sim_timeset", "Seq-Sim"),
        ("deepgo-se", "DeepGO-SE"),
        ("sprof-go", "SPROF-GO"),
    ]:
        p = REV2_DIR / f"stage2_and_baselines_{m_key}.json"
        if not p.exists():
            continue
        cc = load_json(p).get("subontologies", {}).get("cc", {})
        for s_key, s_label in spec_order:
            if s_key not in cc:
                continue
            r = cc[s_key]
            g50 = r.get("gaef_all50", {})
            g41 = r.get("gaef_gt_subset", {})
            cafa = r.get("pooled_cafa", {})
            lines.append(
                f"| {m_label} | {s_label} | {g50.get('taxon_consistency', 0.0):.2f} | "
                f"{100.0 * g50.get('complex_coherence', 0.0):.2f} | "
                f"{g41.get('taxon_consistency', 0.0):.2f} | "
                f"{100.0 * g41.get('complex_coherence', 0.0):.2f} | "
                f"{cafa.get('fmax', 0.0):.4f} | {cafa.get('smin', 0.0):.2f} | "
                f"{cafa.get('avg_auc', 0.0):.4f} | {cafa.get('aupr', 0.0):.4f} | "
                f"{cafa.get('wfmax', 0.0):.4f} |"
            )
    return "\n".join(lines) + "\n"


def build_table_s15_complex_portal() -> str:
    """Reproduce Supplementary Table S15 (Complex Portal validation of Standard vs Demotion-Only Stage 2)."""
    lines = [
        "# Supplementary Table S15: Complex Portal Subunit Validation of Stage 2 Modes",
        "",
        "| Predictor | Scope (n) | Stage 2 Mode | Total Prom. | Covered Prom. | True Subunit | Total Dem. | Covered Dem. | Non-Subunit Removal |",
        "| :--- | :--- | :--- | ---: | ---: | :--- | ---: | ---: | :--- |",
    ]
    configs = [
        ("DeepGO-SE", "Timeset (50)", "Standard (taxons)", "cp_val_deepgo-se_standard_withtax.json"),
        ("DeepGO-SE", "Timeset (50)", "Demotion-Only (taxons)", "cp_val_deepgo-se_demotion_only_withtax.json"),
        ("DeepGO-SE", "Timeset (50)", "Demotion-Only (no taxons)", "cp_val_deepgo-se_demotion_only_notax.json"),
        ("DeepGO-SE", "CP Ref. (4)", "Demotion-Only (no taxons)", "cp6_val_deepgo-se_demotion_only.json"),
        ("Seq-Sim", "Timeset (50)", "Standard (taxons)", "cp_val_seq_sim_timeset_standard_withtax.json"),
        ("Seq-Sim", "Timeset (50)", "Demotion-Only (taxons)", "cp_val_seq_sim_timeset_demotion_only_withtax.json"),
        ("Seq-Sim", "Timeset (50)", "Demotion-Only (no taxons)", "cp_val_seq_sim_timeset_demotion_only_notax.json"),
        ("SPROF-GO", "Timeset (50)", "Standard (taxons)", "cp_val_sprof-go_standard_withtax.json"),
        ("SPROF-GO", "Timeset (50)", "Demotion-Only (taxons)", "cp_val_sprof-go_demotion_only_withtax.json"),
        ("SPROF-GO", "Timeset (50)", "Demotion-Only (no taxons)", "cp_val_sprof-go_demotion_only_notax.json"),
    ]
    for pred, scope, mode, fname in configs:
        p = REV2_DIR / fname
        if not p.exists():
            continue
        d = load_json(p)
        pr = d["promotions"]
        dm = d["demotions"]
        p_cov = pr["with_curated_reference"]
        p_true = pr["true_subunit_promotions"]
        p_str = f"{p_true} / {p_cov} ({pr['precision_pct_on_covered']:.1f}%)" if p_cov > 0 else "--"
        d_cov = dm["with_curated_reference"]
        d_non = d_cov - dm["harmful_demotions_curated_subunit_removed"]
        d_str = f"{d_non} / {d_cov} ({dm['good_demotion_pct_on_covered']:.1f}%)" if d_cov > 0 else "--"
        lines.append(
            f"| {pred} | {scope} | {mode} | {pr['total']:,} | {p_cov} | {p_str} | "
            f"{dm['total']:,} | {d_cov} | {d_str} |"
        )
    return "\n".join(lines) + "\n"


def copy_other_revision_artifacts() -> list[Path]:
    """Copy committed TSV/CSV/JSON summary tables and figures into output/."""
    written = []
    paper_fig_dir = REPO_ROOT.parent / "paper" / "oup-authoring-template" / "figures"
    mapping = {
        "TableS6_ablation_var_counts.csv": REV_DIR / "ablation_var_counts.csv",
        "TableS7_top20_violated_pairs.tsv": REV_DIR / "top20_violated_pairs.tsv",
        "TableS11_cost_ablation.txt": REV_DIR / "cost_ablation.txt",
        "TableS14_soft_regularizer_deepgo-se_cc.tsv": REV_DIR / "soft_regularizer_deepgo-se_cc.tsv",
        "FigureS1_metagenomics_bp.pdf": paper_fig_dir / "metagenomics_fig_main_condense.pdf",
        "FigureS2_metagenomics_seqsim.pdf": paper_fig_dir / "metagenomics_fig_summary.pdf",
    }
    for dst_name, src_path in mapping.items():
        if src_path.exists():
            dst = OUT_DIR / dst_name
            shutil.copy2(src_path, dst)
            written.append(dst)
    return written


def extract_all_latex_tables() -> list[Path]:
    """Extract all LaTeX tables from paper-draft.tex (Tables 1-2) and paper-supplementary.tex (Tables S1-S15) if present."""
    written = []
    paper_dir = REPO_ROOT.parent / "paper" / "oup-authoring-template"
    for tex_name, prefix in [("paper-draft.tex", "Main_Table"), ("paper-supplementary.tex", "Supp_Table_S")]:
        tex_path = paper_dir / tex_name
        if not tex_path.exists():
            continue
        text = tex_path.read_text(encoding="utf-8")
        blocks = re.findall(r"(\\begin\{table\*?\}.*?\\end\{table\*?\})", text, flags=re.DOTALL)
        for idx, block in enumerate(blocks, start=1):
            out = write_artifact(f"{prefix}{idx}.tex", block + "\n")
            written.append(out)
    return written


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated = [
        write_artifact("Table1_flips.md", build_table1_flips()),
        write_artifact("TableS10_tau_sweep.md", build_table_s10_tau_sweep()),
        write_artifact("TableS12_lineage_inference.md", build_table_s12_lineage()),
        write_artifact("TableS13_baselines_soft_stage2.md", build_table_s13_baselines_and_stage2()),
        write_artifact("TableS15_complex_portal.md", build_table_s15_complex_portal()),
    ]
    generated.extend(copy_other_revision_artifacts())
    generated.extend(extract_all_latex_tables())

    print(f"Successfully generated {len(generated)} reproducibility artifacts in {OUT_DIR}:")
    for p in generated:
        print(f"  - {p.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
