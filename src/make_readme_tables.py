"""Print the README tables from the CSVs in results/ — run from the repo root.

    python src/make_readme_tables.py

The README quotes numbers only through these tables, so README and results/ cannot drift
apart again (which is exactly what happened in the first iteration).
"""
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1] / "results"


def table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(out)


def f3(x):
    return "—" if pd.isna(x) else f"{x:.3f}"


def main():
    # ---- Table 1: single-run all-cell metrics (notebook 05) ---------------------------
    best = pd.read_csv(R / "best_per_method.csv")
    t1 = pd.DataFrame({"Method": best["method"], "Leiden res.": best["resolution"],
                       "ARI (all cells)": best["ARI"].map(f3), "NMI (all cells)": best["NMI"].map(f3),
                       "clusters": best["n_clusters"]})
    print("### Table 1 — single run, all cells (notebook 05)\n")
    print(table(t1), "\n")

    # ---- Table 2: five seeds (notebook 07) ---------------------------------------------
    s = pd.read_csv(R / "seed_stability.csv")
    rows = []
    for m in ["scVI", "scANVI"]:
        d = s[s["method"] == m]
        row = {"Method": m, "seeds": len(d),
               "ARI all cells": f"{d['ARI_all'].mean():.3f} ± {d['ARI_all'].std():.3f}",
               "NMI all cells": f"{d['NMI_all'].mean():.3f} ± {d['NMI_all'].std():.3f}"}
        if m == "scANVI":
            row["held-out accuracy"] = f"{d['accuracy_heldout'].mean():.3f} ± {d['accuracy_heldout'].std():.3f}"
            row["held-out macro-F1"] = f"{d['macroF1_heldout'].mean():.3f} ± {d['macroF1_heldout'].std():.3f}"
        else:
            row["held-out accuracy"] = "—"; row["held-out macro-F1"] = "—"
        rows.append(row)
    print("### Table 2 — five seeds, mean ± sd (notebook 07)\n")
    print(table(pd.DataFrame(rows)), "\n")
    ep = s.loc[s["method"] == "scVI", "epochs"].tolist()
    print(f"scVI epochs trained per seed: {ep} (early stopping fired: {[e < 200 for e in ep]})\n")

    # ---- Table 3: held-out views (notebook 05 §4b) --------------------------------------
    h = pd.read_csv(R / "heldout_metrics.csv")
    t3 = pd.DataFrame({"Method": h["method"],
                       "ARI all cells": h["ARI_all"].map(f3),
                       "ARI held-out only": h["ARI_heldout"].map(f3),
                       "kNN transfer acc. → smartseq2": h["kNN_transfer_acc"].map(f3),
                       "kNN transfer macro-F1": h["kNN_transfer_macroF1"].map(f3),
                       "batch mixing (50-NN)": h["batch_mixing_all"].map(f3),
                       "cell-type purity (50-NN)": h["celltype_purity_all"].map(f3)})
    print("### Table 3 — held-out views and integration axes (notebook 05 §4b)\n")
    print(table(t3), "\n")

    # ---- Table 4: error anatomy verdict inputs (notebook 08) ---------------------------
    p = R / "error_anatomy_pairs.csv"
    if p.exists():
        e = pd.read_csv(p)
        t4 = pd.DataFrame({
            "Error (true → predicted)": e["pair"],
            "stable / unstable cells": [f"{a} / {b}" for a, b in zip(e["n_stable"], e["n_unstable"])],
            "errors with predicted-type marker above p95 of true type": e["errors with pred_marker > p95 of true type"],
            "errors with true-type marker below p05 of true type": e["errors with true_marker < p05 of true type"],
            "doublet score, errors vs correct (median)": [f"{a:.3f} vs {b:.3f}" for a, b in
                                                          zip(e["doublet score, median errors"], e["doublet score, median correct true type"])],
            "genes detected, errors vs correct (median)": [f"{int(a)} vs {int(b)}" for a, b in
                                                            zip(e["n_genes, median errors"], e["n_genes, median correct true type"])],
            "mean confidence": e["mean confidence (stable errors)"].map(lambda x: f"{x:.2f}"),
        })
        print("### Table 4 — stable label-transfer errors on smartseq2, evidence per pair (notebook 08)\n")
        print(table(t4), "\n")

    m = R / "error_anatomy_harmony_merges.csv"
    if m.exists():
        g = pd.read_csv(m)
        t5 = pd.DataFrame({"Cells": [f"{a} in the {b} cluster" for a, b in zip(g["celltype"], g["cluster_dominant"])],
                           "n": g["n"], "top technology": g["top_tech"],
                           "share from top technology": [f"{a:.0%} (expected {b:.0%})" for a, b in
                                                         zip(g["share_top_tech"], g["share_top_tech_expected"])]})
        print("### Table 5 — Harmony merges at resolution 0.3: where the 'wrong' cells come from (notebook 08, Part A)\n")
        print(table(t5), "\n")


if __name__ == "__main__":
    main()
