import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT_FILE = r"C:\Users\91892\kalralabintern\residue_count\protein_with_aromatic_residue_counts.csv"

OUTPUT_DIR = r"C:\Users\91892\kalralabintern\residue_count\Aromatic_Density_Analysis"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE)

# ============================================================
# CALCULATE DENSITIES
# ============================================================

df["W_density"] = df["Tryptophan(W)"] / df["Length"]
df["Y_density"] = df["Tyrosine(Y)"] / df["Length"]
df["F_density"] = df["Phenylalanine(F)"] / df["Length"]

df["Total_Aromatic"] = (
    df["Tryptophan(W)"] +
    df["Tyrosine(Y)"] +
    df["Phenylalanine(F)"]
)

df["Total_Aromatic_Density"] = (
    df["Total_Aromatic"] /
    df["Length"]
)

# Save updated table
df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "protein_with_aromatic_densities.csv"
    ),
    index=False
)

# ============================================================
# SUMMARY STATISTICS
# ============================================================

summary = df[
    [
        "W_density",
        "Y_density",
        "F_density",
        "Total_Aromatic_Density"
    ]
].describe()

summary.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "Density_summary_statistics.csv"
    )
)

print(summary)

# ============================================================
# VIOLIN PLOT
# ============================================================

plt.figure(figsize=(7,6))

data = [
    df["W_density"],
    df["Y_density"],
    df["F_density"],
    df["Total_Aromatic_Density"]
]

labels = [
    "W",
    "Y",
    "F",
    "Total"
]

plt.violinplot(
    data,
    showmeans=True,
    showmedians=True
)

plt.xticks(
    np.arange(1,5),
    labels,
    fontsize=12
)

plt.ylabel("Residues per amino acid")
plt.title("Aromatic Residue Density")

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Violin_density.png"
    ),
    dpi=600
)

plt.close()

# ============================================================
# BOXPLOT
# ============================================================

plt.figure(figsize=(7,6))

plt.boxplot(
    data,
    tick_labels=labels
)

plt.ylabel("Residues per amino acid")
plt.title("Aromatic Residue Density")

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Boxplot_density.png"
    ),
    dpi=600
)

plt.close()

# ============================================================
# HISTOGRAMS
# ============================================================

cols = [
    "W_density",
    "Y_density",
    "F_density",
    "Total_Aromatic_Density"
]

for col in cols:

    plt.figure(figsize=(6,5))

    plt.hist(
        df[col],
        bins=30
    )

    plt.xlabel(col)
    plt.ylabel("Number of Proteins")
    plt.title(col)

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"{col}_histogram.png"
        ),
        dpi=600
    )

    plt.close()

# ============================================================
# CORRELATION HEATMAP
# ============================================================

corr = df[
    [
        "Length",
        "W_density",
        "Y_density",
        "F_density",
        "Total_Aromatic_Density"
    ]
].corr()

fig, ax = plt.subplots(figsize=(7,6))

im = ax.imshow(corr)

ax.set_xticks(np.arange(len(corr.columns)))
ax.set_xticklabels(corr.columns, rotation=45, ha="right")

ax.set_yticks(np.arange(len(corr.columns)))
ax.set_yticklabels(corr.columns)

for i in range(len(corr)):
    for j in range(len(corr)):
        ax.text(
            j,
            i,
            f"{corr.iloc[i,j]:.2f}",
            ha="center",
            va="center",
            color="white" if abs(corr.iloc[i,j]) > 0.5 else "black"
        )

plt.colorbar(im)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Correlation_heatmap.png"
    ),
    dpi=600
)

plt.close()

print("\nAnalysis completed successfully.")
print(f"Results saved to:\n{OUTPUT_DIR}")