import pandas as pd
import re

# ============================================================
# FILES
# ============================================================

PANTHER_FILE = r"C:\Users\91892\kalralabintern\PTHR19.0_human.txt"

INPUT_CSV = r"C:\Users\91892\kalralabintern\pr_fams\main\diff_vs_protein_families_with_superfamilies.csv"

OUTPUT_CSV = r"C:\Users\91892\kalralabintern\pr_fams\main\diff_vs_protein_families_with_protein_class_panther.csv"

# ============================================================
# READ PANTHER
# ============================================================

panther = pd.read_csv(
    PANTHER_FILE,
    sep="\t",
    header=None,
    low_memory=False
)

panther.columns = [
    "Gene_ID",
    "UniProt_ID",
    "Gene_Symbol",
    "PANTHER_Family",
    "Protein_Name",
    "Description",
    "GO_Molecular_Function",
    "GO_Biological_Process",
    "GO_Cellular_Component",
    "Protein_Class",
    "Pathway"
]

print("Proteins in PANTHER:", len(panther))

# ============================================================
# CLEAN PROTEIN CLASS
# ============================================================

def clean_class(x):

    if pd.isna(x):
        return None

    # membrane traffic protein#PC00150
    x = str(x).split("#")[0]

    return x.strip()

panther["Protein_Class"] = (
    panther["Protein_Class"]
    .apply(clean_class)
)

# ============================================================
# BUILD LOOKUP
# ============================================================

lookup = dict(
    zip(
        panther["UniProt_ID"],
        panther["Protein_Class"]
    )
)

print("Protein classes:", len(lookup))

# ============================================================
# READ YOUR CSV
# ============================================================

df = pd.read_csv(INPUT_CSV)

# ============================================================
# EXTRACT UNIPROT ID
# ============================================================

def extract_uniprot(x):

    if pd.isna(x):
        return None

    x = str(x)

    # AF-Q9Y6K9-F1
    m = re.search(r'AF-([A-Z0-9]+)-F', x)

    if m:
        return m.group(1)

    return x

df["UniProt_ID"] = (
    df["Alphafold_ID"]
    .apply(extract_uniprot)
)

# ============================================================
# MAP PROTEIN CLASS
# ============================================================

df["Protein_Class"] = (
    df["UniProt_ID"]
    .map(lookup)
)

# ============================================================
# REPORT
# ============================================================

mapped = df["Protein_Class"].notna().sum()

print(f"\nMapped {mapped:,} proteins")

print(df["Protein_Class"].value_counts().head(20))

# ============================================================
# SAVE
# ============================================================

df.to_csv(
    OUTPUT_CSV,
    index=False
)

print("\nSaved to:")
print(OUTPUT_CSV)