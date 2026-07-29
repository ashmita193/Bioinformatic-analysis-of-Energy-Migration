import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ==========================================================
# INPUT / OUTPUT
# ==========================================================

INPUT_FILE = r"C:\Users\91892\kalralabintern\residue_count\protein_with_aromatic_residue_counts.csv"

OUTPUT_DIR = r"C:\Users\91892\kalralabintern\residue_count\PCA_Output"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================================
# LOAD DATA
# ==========================================================

df = pd.read_csv(INPUT_FILE)

# ==========================================================
# FEATURES FOR PCA
# ==========================================================

features = [
    "Length",
    "Tryptophan(W)",
    "Tyrosine(Y)",
    "Phenylalanine(F)"
]

X = df[features].copy()

# Remove missing values
X = X.dropna()

# Keep corresponding rows
df = df.loc[X.index].reset_index(drop=True)

# ==========================================================
# STANDARDIZE
# ==========================================================

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

# ==========================================================
# PCA
# ==========================================================

pca = PCA()

scores = pca.fit_transform(X_scaled)

explained = pca.explained_variance_ratio_

# ==========================================================
# SAVE PCA SCORES
# ==========================================================

score_df = pd.DataFrame(
    scores,
    columns=[f"PC{i+1}" for i in range(len(features))]
)

score_df.to_csv(
    os.path.join(OUTPUT_DIR, "PCA_scores.csv"),
    index=False
)

# ==========================================================
# SAVE LOADINGS
# ==========================================================

loadings = pd.DataFrame(
    pca.components_.T,
    index=features,
    columns=[f"PC{i+1}" for i in range(len(features))]
)

loadings.to_csv(
    os.path.join(OUTPUT_DIR, "PCA_loadings.csv")
)

# ==========================================================
# SAVE EXPLAINED VARIANCE
# ==========================================================

variance_df = pd.DataFrame({
    "Principal Component":
        [f"PC{i+1}" for i in range(len(features))],

    "Explained Variance":
        explained,

    "Explained Variance (%)":
        explained*100
})

variance_df.to_csv(
    os.path.join(OUTPUT_DIR, "Explained_variance.csv"),
    index=False
)

# ==========================================================
# PRINT SUMMARY
# ==========================================================

print("\nExplained Variance\n")

for i, v in enumerate(explained):
    print(f"PC{i+1}: {v*100:.2f}%")

print(f"\nPC1 + PC2 = {(explained[:2].sum()*100):.2f}%")

print("\nLoadings\n")
print(loadings)

# ==========================================================
# SCREE PLOT
# ==========================================================

plt.figure(figsize=(7,5))

plt.plot(
    range(1, len(features)+1),
    explained*100,
    marker='o',
    linewidth=2
)

plt.xticks(range(1, len(features)+1))

plt.xlabel("Principal Component", fontsize=13)
plt.ylabel("Variance Explained (%)", fontsize=13)
plt.title("Scree Plot", fontsize=16)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "Scree_plot.png"),
    dpi=600
)

plt.close()

# ==========================================================
# PCA SCORE PLOT
# ==========================================================

plt.figure(figsize=(7,6))

plt.scatter(
    score_df["PC1"],
    score_df["PC2"],
    s=35,
    alpha=0.75
)

plt.xlabel(
    f"PC1 ({explained[0]*100:.1f}%)",
    fontsize=13
)

plt.ylabel(
    f"PC2 ({explained[1]*100:.1f}%)",
    fontsize=13
)

plt.title("PCA Score Plot", fontsize=16)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "PCA_score_plot.png"),
    dpi=600
)

plt.close()

# ==========================================================
# LOADING PLOT
# ==========================================================

plt.figure(figsize=(8,8))

for feature in features:

    x = loadings.loc[feature, "PC1"]
    y = loadings.loc[feature, "PC2"]

    plt.arrow(
        0,
        0,
        x,
        y,
        head_width=0.03,
        color="tab:blue",
        length_includes_head=True
    )

    plt.text(
        x*1.08,
        y*1.08,
        feature,
        fontsize=12
    )

plt.axhline(0, color='gray')
plt.axvline(0, color='gray')

plt.xlim(-1,1)
plt.ylim(-1,1)

plt.xlabel("PC1", fontsize=13)
plt.ylabel("PC2", fontsize=13)

plt.title("PCA Loading Plot", fontsize=16)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "PCA_loading_plot.png"),
    dpi=600
)

plt.close()

# ==========================================================
# PCA BIPLOT
# ==========================================================

plt.figure(figsize=(8,8))

# Scores
plt.scatter(
    score_df["PC1"],
    score_df["PC2"],
    alpha=0.35,
    s=20
)

# Scale arrows
scale = 3

for feature in features:

    x = loadings.loc[feature, "PC1"] * scale
    y = loadings.loc[feature, "PC2"] * scale

    plt.arrow(
        0,
        0,
        x,
        y,
        color="red",
        head_width=0.08,
        length_includes_head=True
    )

    plt.text(
        x*1.05,
        y*1.05,
        feature,
        fontsize=12,
        color="red"
    )

plt.axhline(0, color="gray")
plt.axvline(0, color="gray")

plt.xlabel(
    f"PC1 ({explained[0]*100:.1f}%)",
    fontsize=13
)

plt.ylabel(
    f"PC2 ({explained[1]*100:.1f}%)",
    fontsize=13
)

plt.title("PCA Biplot", fontsize=16)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "PCA_biplot.png"),
    dpi=600
)

plt.close()

print("\n========================================")
print("PCA analysis completed successfully.")
print("Results saved to:")
print(OUTPUT_DIR)
print("========================================")