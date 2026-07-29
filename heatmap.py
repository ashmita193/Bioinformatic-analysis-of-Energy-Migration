import pandas as pd
import numpy as np
import scipy.stats as stats
import scikit_posthocs as sp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# ==========================================================
# INPUT
# ==========================================================

INPUT_FILE = r"C:\Users\91892\kalralabintern\cleaned_det\diff_det_fret_vs_prfams_broadclasses.csv"

GROUP_COLUMN = "Final_Broad_Class"
VALUE_COLUMN = "Diffusion_Length_DET"

P_ADJUST = "holm"      # bonferroni, holm, fdr_bh etc.
ALPHA = 0.05

# ==========================================================
# LOAD DATA
# ==========================================================

df = pd.read_csv(INPUT_FILE)

df = df[[GROUP_COLUMN, VALUE_COLUMN]].dropna()

# ==========================================================
# KRUSKAL-WALLIS TEST
# ==========================================================

groups = [
    group[VALUE_COLUMN].values
    for _, group in df.groupby(GROUP_COLUMN)
]

H, p = stats.kruskal(*groups)

print("=" * 60)
print("KRUSKAL-WALLIS TEST")
print("=" * 60)
print(f"H statistic : {H:.4f}")
print(f"P-value     : {p:.3e}")

if p < ALPHA:
    print("Result      : Significant")
else:
    print("Result      : Not Significant")

# ==========================================================
# DUNN POST-HOC TEST
# ==========================================================

dunn = sp.posthoc_dunn(
    df,
    val_col=VALUE_COLUMN,
    group_col=GROUP_COLUMN,
    p_adjust=P_ADJUST
)

dunn.to_csv("Dunn_Posthoc_Pvalues.csv")

print("\nPairwise adjusted p-values saved.")

# ==========================================================
# PREPARE HEATMAP
# ==========================================================

labels = dunn.columns

n = len(labels)

# lower triangle only
mask = np.triu(np.ones((n, n), dtype=bool))

# significance matrix
sig_matrix = (dunn < ALPHA).astype(int)

# hide upper triangle
sig_matrix = sig_matrix.mask(mask)

# annotation text
annot = dunn.copy().astype(object)

for i in range(n):
    for j in range(n):

        if mask[i, j]:
            annot.iloc[i, j] = ""

        elif i == j:
            annot.iloc[i, j] = "—"

        else:

            pval = dunn.iloc[i, j]

            if pval < 1e-3:
                annot.iloc[i, j] = f"{pval:.0e}"
            else:
                annot.iloc[i, j] = f"{pval:.3f}"

# ==========================================================
# PLOT
# ==========================================================

cmap = ListedColormap([
    "#d9d9d9",    # Not significant
    "#2ca25f"     # Significant
])

fig, ax = plt.subplots(figsize=(9,8))

im = ax.imshow(sig_matrix, cmap=cmap, vmin=0, vmax=1)

# ticks
ax.set_xticks(np.arange(n))
ax.set_yticks(np.arange(n))

ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=11)
ax.set_yticklabels(labels, fontsize=11)

# grid
ax.set_xticks(np.arange(-.5, n, 1), minor=True)
ax.set_yticks(np.arange(-.5, n, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=2)

ax.tick_params(which="minor", bottom=False, left=False)

# annotations
for i in range(n):
    for j in range(n):

        if annot.iloc[i, j] != "":

            color = "white" if sig_matrix.iloc[i, j] == 1 else "black"

            ax.text(
                j,
                i,
                annot.iloc[i, j],
                ha="center",
                va="center",
                fontsize=10,
                color=color,
                fontweight="bold"
            )

# legend
legend = [
    Patch(facecolor="#2ca25f", edgecolor="none",
          label=f"Significant (p < {ALPHA})"),
    Patch(facecolor="#d9d9d9", edgecolor="none",
          label=f"Not significant (p ≥ {ALPHA})")
]

ax.legend(handles=legend,
          loc="upper right",
          frameon=True)

'''plt.title(
    f"Dunn's Post Hoc Test ({P_ADJUST.capitalize()} correction)",
    fontsize=15,
    weight="bold"
)'''


plt.savefig(
    r"C:\Users\91892\kalralabintern\cleaned_det\Dunn_Posthoc_Heatmap_DET.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()