"""Standalone runner for the new section 4b of 05_evaluation.ipynb.

Everything above the marker reproduces the state 05 has after its section 4
(adata_baseline with X_pca and leiden_r* columns, method_files, RESOLUTIONS, best_per_method).
Everything between the markers is the cell to paste into 05 verbatim.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sc.settings.verbosity = 1
RESOLUTIONS = [0.3, 0.5, 0.7, 1.0, 1.3]

adata_baseline = sc.read("../data/raw/pancreas.h5ad")
sc.pp.highly_variable_genes(adata_baseline, n_top_genes=2000, batch_key="tech", subset=False)
adata_baseline = adata_baseline[:, adata_baseline.var["highly_variable"]].copy()
sc.pp.scale(adata_baseline, max_value=10)
sc.tl.pca(adata_baseline, n_comps=30)
sc.pp.neighbors(adata_baseline, use_rep="X_pca")
for res in RESOLUTIONS:
    sc.tl.leiden(adata_baseline, resolution=res, flavor="igraph", n_iterations=2, directed=False,
                 key_added=f"leiden_r{res}")

method_files = {
    "Harmony": "../data/processed/pancreas_harmony.h5ad",
    "scVI": "../data/processed/pancreas_scvi.h5ad",
    "scANVI": "../data/processed/pancreas_scanvi.h5ad",
}
best_per_method = pd.read_csv("../results/best_per_method.csv")

# ======================= BEGIN cell to insert into 05, section 4b ========================
# The table above has a flaw: scANVI was trained with the labels of ~85 % of the cells, so
# its all-cell ARI/NMI is not comparable with Harmony or scVI. The natural fix — score only
# the 2394 held-out smartseq2 cells — turns out to be uninformative for a different reason,
# so this section reports three views and says what each one can and cannot show.
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.metrics import accuracy_score, f1_score

HOLDOUT_TECH = "smartseq2"
latents = {"baseline (no integration)": (adata_baseline, "X_pca")}
for m, p in method_files.items():
    ad_m = sc.read(p)
    latents[m] = (ad_m, {"Harmony": "X_pca_harmony", "scVI": "X_scVI", "scANVI": "X_scANVI"}[m])

K = 50
rows = []
for method_name, (ad_m, rep) in latents.items():
    obs, Z = ad_m.obs, ad_m.obsm[rep]
    tech = obs["tech"].astype(str).values
    ct = obs["celltype"].astype(str).values
    ho = tech == HOLDOUT_TECH
    res = float(best_per_method.loc[best_per_method["method"] == method_name, "resolution"].iloc[0])
    cl = obs[f"leiden_r{res}"]

    # (1) Leiden ARI/NMI restricted to the held-out cells
    ari_ho = adjusted_rand_score(ct[ho], cl[ho])
    nmi_ho = normalized_mutual_info_score(ct[ho], cl[ho])

    # (2) cross-batch kNN label transfer: reference = all labelled non-smartseq2 cells
    knn = KNeighborsClassifier(n_neighbors=15, weights="distance").fit(Z[~ho], ct[~ho])
    pred = knn.predict(Z[ho])
    acc = accuracy_score(ct[ho], pred)
    f1 = f1_score(ct[ho], pred, average="macro", labels=sorted(set(ct[ho])))  # same label set as in 07

    # (3) batch mixing vs cell-type purity of the K nearest neighbours (all cells)
    idx = NearestNeighbors(n_neighbors=K + 1).fit(Z).kneighbors(Z, return_distance=False)[:, 1:]
    other_tech = (tech[idx] != tech[:, None]).mean(1)
    expected = 1 - pd.Series(tech).value_counts(normalize=True)[tech].values
    mixing = np.minimum(other_tech / expected, 1)           # 1 = neighbours as mixed as random
    purity = (ct[idx] == ct[:, None]).mean(1)               # 1 = neighbours all same cell type

    rows.append({"method": method_name, "resolution": res,
                 "ARI_heldout": ari_ho, "NMI_heldout": nmi_ho,
                 "kNN_transfer_acc": acc, "kNN_transfer_macroF1": f1,
                 "batch_mixing_all": mixing.mean(), "batch_mixing_heldout": mixing[ho].mean(),
                 "celltype_purity_all": purity.mean(), "celltype_purity_heldout": purity[ho].mean()})

heldout = (best_per_method[["method", "ARI", "NMI"]].rename(columns={"ARI": "ARI_all", "NMI": "NMI_all"})
           .merge(pd.DataFrame(rows), on="method"))
heldout.to_csv("../results/heldout_metrics.csv", index=False)
print(heldout.round(3).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
x = np.arange(len(heldout)); w = 0.27
axes[0].bar(x - w, heldout["ARI_all"], w, label="ARI, all cells (table above)")
axes[0].bar(x, heldout["ARI_heldout"], w, label="ARI, held-out smartseq2 only")
axes[0].bar(x + w, heldout["kNN_transfer_acc"], w, label="kNN label transfer → smartseq2")
axes[0].set_xticks(x); axes[0].set_xticklabels(heldout["method"], rotation=12, fontsize=8)
axes[0].set_ylim(0, 1.05); axes[0].legend(fontsize=8); axes[0].set_title("Held-out views: all near the ceiling")
for _, r in heldout.iterrows():
    axes[1].scatter(r["batch_mixing_all"], r["celltype_purity_all"], s=80)
    axes[1].annotate(r["method"], (r["batch_mixing_all"], r["celltype_purity_all"]),
                     textcoords="offset points", xytext=(6, 4), fontsize=8)
axes[1].set_xlabel(f"batch mixing of {K}-NN (1 = as mixed as random)")
axes[1].set_ylabel(f"cell-type purity of {K}-NN")
axes[1].set_xlim(0, 1.05); axes[1].set_ylim(0.90, 1.0); axes[1].set_title("What integration actually changes")
axes[1].grid(alpha=0.3)
plt.tight_layout()
plt.savefig("../results/figures/05_heldout_views.png", dpi=150)
plt.show()
# ======================== END cell to insert into 05, section 4b =========================
