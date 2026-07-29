import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import textwrap
from matplotlib.patches import Patch
# ============================================================
# INPUT
# ============================================================

INPUT_FILE = r"C:\Users\91892\kalralabintern\cleaned_det\diff_det_fret_vs_prfams_broadclasses.csv"
df = pd.read_csv(INPUT_FILE)

# ============================================================
# CLEAN DATA
# ============================================================

df = df.dropna(subset=["Final_Broad_Class", "Diffusion_Length_DET"])

df["Diffusion_Length_DET"] = pd.to_numeric(
    df["Diffusion_Length_DET"],
    errors="coerce"
)

df = df.dropna(subset=["Diffusion_Length_DET"])

# ============================================================
# STATISTICS
# ============================================================

# ============================================================
# STATISTICS + CUSTOM ORDER
# ============================================================

stats = (
    df.groupby("Final_Broad_Class")
      .agg(
          Mean=("Diffusion_Length_DET", "mean"),
          SD=("Diffusion_Length_DET", "std"),
          N=("Diffusion_Length_DET", "count"),
          n=("Alphafold_ID", "nunique")
      )
)

# Define the exact left-to-right order you want
order = [
   "Transcription Regulation",
   "Immune",
   "Membrane Traffic Protein",
   "Structural",
   "Cell Adhesion",
   "Chaperone",
   "Nucleic Acid & Translation",
   "Signalling",
   "Transporter / Channel",
   "Receptor",
   "Viral or Transposable Element",
   "Enzyme"
    # Add all your remaining classes here
]

# Keep only classes that are actually present
order = [cls for cls in order if cls in stats.index]

# Reorder statistics
stats = stats.loc[order]

plot_data = [
    df.loc[
        df["Final_Broad_Class"] == c,
        "Diffusion_Length_DET"
    ].values
    for c in order
]

# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(figsize=(22, 9))

positions = np.arange(len(order)) * 1.8 + 1
width = 0.76

class_colors = {
    "Transcription Regulation": "#E5C8A8",
    "Immune": "#D9E5F5",
    "Membrane Traffic Protein": "#FFD9A6",
    "Structural": "#C9DEF7",
    "Cell Adhesion": "#D8EFBF",
    "Chaperone": "#F7D8E8",
   "Nucleic Acid & Translation": "#ECECEC",
    "Signalling": "#DDC9F0",
   "Transporter / Channel": "#D9EFD4",
   "Receptor": "#FBEAA7",
   "Viral or Transposable Element": "#CBD7EE" ,
   "Enzyme": "#F4C2C2"
    
    # Add the remaining classes and colors
}

# ============================================================
# DRAW VIOLINS
# ============================================================

for pos, vals, color in zip(
        positions,
        plot_data,
        [class_colors.get(cls) for cls in order]):

    vals = np.asarray(vals)

    if len(vals) < 5:
        continue

    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1

    lower_limit = q1 - 1.5 * iqr
    upper_limit = q3 + 1.5 * iqr

    low = vals[vals >= lower_limit].min()
    high = vals[vals <= upper_limit].max()

    kde = gaussian_kde(vals, bw_method=0.30)

    y = np.linspace(vals.min(), vals.max(), 600)

    density = kde(y)
    density = density / density.max() * width

    mask = (y >= low) & (y <= high)

    y = y[mask]
    density = density[mask]

    ax.fill_betweenx(
        y,
        pos - density,
        pos + density,
        facecolor=color,
        edgecolor="#444444",
        linewidth=1.5,
        alpha=0.75,
        zorder=1
    )

# ============================================================
# BOXPLOTS
# ============================================================

ax.boxplot(
    plot_data,
    positions=positions,
    widths=0.60,
    patch_artist=True,
    showfliers=False,
    whis=1.5,
    boxprops=dict(
    facecolor="white",
    edgecolor="black",
    linewidth=1.5
    ),
    whiskerprops=dict(color="black"),
    capprops=dict(color="black"),
    medianprops=dict(
        color="#0072B2",
        linewidth=2.3
    )
)

# ============================================================
# MEANS
# ============================================================

ax.scatter(
    positions,
    stats["Mean"],
    marker="D",
    color="#D62728",
    edgecolor="darkred",
    s=60,
    zorder=5
)

# ============================================================
# LABELS
# ============================================================

labels = []

for cls in order:

    row = stats.loc[cls]

    wrapped = textwrap.fill(cls, width=18)

    labels.append(
        f"{wrapped}\n"
        f"N={int(row['N'])}\n"
        f"μ={row['Mean']:.2f}\n"
        f"σ={row['SD']:.2f}"
    )

ax.set_xticks(positions)
ax.set_xticklabels(
    labels,
    fontsize=14,
    fontfamily="Sans-serif"
)

ax.set_xlim(
    positions[0] - 1,
    positions[-1] + 1
)

ax.set_xlabel(
    "Broad Protein Class",
    fontsize=18,
    fontweight="bold",
    fontfamily="Arial"
)

ax.set_ylabel(
    "Diffusion Length (Å)",
    fontsize=18,
    fontweight="bold",
    fontfamily="Arial"
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=1.5,
    alpha=0.35
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ============================================================
# LEGEND
# ============================================================

legend = [
    Patch(
        facecolor="white",
        edgecolor="black",
        label="IQR (25th–75th percentile)"
    ),
    Line2D(
        [0],
        [0],
        color="#0072B2",
        lw=2.3,
        label="Median"
    ),
    Line2D(
        [0],
        [0],
        color="black",
        lw=1,
        label="Whiskers (1.5×IQR)"
    ),
    Line2D(
        [0],
        [0],
        marker="D",
        color="w",
        markerfacecolor="#D62728",
        markeredgecolor="darkred",
        markersize=8,
        label="Mean"
    )
]

ax.legend(
    handles=legend,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.14),
    ncol=4,
    frameon=True,
    edgecolor="black",
    facecolor="white",
    fontsize=16,
    framealpha=1
)

plt.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.30,
    top=0.84
)

plt.savefig(
   r"C:\Users\91892\kalralabintern\cleaned_det\new_det.png",
    dpi=600,
    bbox_inches="tight"
)


plt.show()