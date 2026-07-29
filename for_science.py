import pandas as pd

# ============================================================
# INPUT FILES
# ============================================================

PROTEIN_FILE = r"C:\Users\91892\kalralabintern\pr_fams\main\panther_protein_class.xlsx"

CLASS_FILE = r"C:\Users\91892\kalralabintern\pr_fams\main\Panther_ontology\Protein_Class_19.0.txt"

REL_FILE = r"C:\Users\91892\kalralabintern\pr_fams\main\Panther_ontology\Protein_class_relationship.txt"

OUTPUT_FILE = r"C:\Users\91892\kalralabintern\pr_fams\main\panther_classes_with_hierarchy.csv"

# ============================================================
# READ PANTHER CLASS FILE
# ============================================================

classes = pd.read_csv(
    CLASS_FILE,
    sep="\t",
    comment="!",
    header=None,
    names=[
        "PC_ID",
        "Hierarchy_Code",
        "Class_Name",
        "Description"
    ],
    usecols=[0,1,2,3],
    engine="python",
    dtype=str
)

id_to_name = dict(zip(classes.PC_ID, classes.Class_Name))
name_to_id = dict(zip(classes.Class_Name.str.strip(), classes.PC_ID))

# ============================================================
# READ RELATIONSHIP FILE
# ============================================================

rels = pd.read_csv(
    REL_FILE,
    sep="\t",
    comment="!",
    header=None,
    names=[
        "Child_ID",
        "Child_Name",
        "Parent_ID",
        "Parent_Name",
        "Depth"
    ],
    usecols=[0,1,2,3,4],
    engine="python",
    dtype=str
)

parent_of = dict(zip(rels.Child_ID, rels.Parent_ID))

# ============================================================
# FUNCTIONS
# ============================================================

def lineage(pcid):

    chain=[]

    while pd.notna(pcid):

        chain.append(pcid)

        if pcid not in parent_of:
            break

        pcid=parent_of[pcid]

    return chain


def parent_name(pcid):

    if pd.isna(pcid):
        return None

    if pcid not in parent_of:
        return None

    return id_to_name[parent_of[pcid]]


def top_level(pcid):

    if pd.isna(pcid):
        return None

    L=lineage(pcid)

    return id_to_name[L[-1]]


def hierarchy_string(pcid):

    if pd.isna(pcid):
        return None

    return " → ".join(
        id_to_name[x]
        for x in lineage(pcid)
    )

# ============================================================
# FINAL 12 CLASS SLIM
# ============================================================

def broad_class(hierarchy):

    if hierarchy is None:
        return "Unknown"

    h=hierarchy.lower()

    if any(x in h for x in [
        "metabolite interconversion enzyme",
        "protein modifying enzyme",
        "transferase",
        "hydrolase",
        "oxidoreductase",
        "lyase",
        "ligase",
        "isomerase",
        "kinase",
        "protease",
        "phosphatase"
    ]):
        return "Enzyme"

    elif any(x in h for x in [
        "transmembrane signal receptor",
        "receptor"
    ]):
        return "Receptor"

    elif any(x in h for x in [
        "transporter",
        "ion channel",
        "carrier"
    ]):
        return "Transporter / Channel"

    elif any(x in h for x in [
        "gene-specific transcriptional regulator",
        "transcription factor",
        "transcription cofactor"
    ]):
        return "Transcription Regulation"

    elif any(x in h for x in [
        "dna metabolism",
        "rna metabolism",
        "translational protein",
        "ribosomal protein",
        "translation factor"
    ]):
        return "Nucleic Acid & Translation"

    elif any(x in h for x in [
        "protein-binding activity modulator",
        "g-protein",
        "g-protein modulator",
        "scaffold/adaptor"
    ]):
        return "Signalling"

    elif any(x in h for x in [
        "cytoskeletal protein",
        "extracellular matrix",
        "structural protein",
        "cell junction"
    ]):
        return "Structural"

    elif any(x in h for x in [
        "cell adhesion"
    ]):
        return "Cell Adhesion"

    elif any(x in h for x in [
        "motor"
    ]):
        return "Motor Protein"

    elif any(x in h for x in [
        "chaperone"
    ]):
        return "Chaperone"

    elif any(x in h for x in [
        "defense/immunity",
        "cytokine",
        "growth factor",
        "intercellular signal molecule",
        "peptide hormone"
    ]):
        return "Immune / Signalling Molecule"

    elif any(x in h for x in [
        "storage protein",
        "transfer/carrier protein"
    ]):
        return "Storage / Carrier"

    else:
        return "Other"

# ============================================================
# LOAD YOUR DATA
# ============================================================

df=pd.read_excel(PROTEIN_FILE)

COLUMN="Protein_Class"      # CHANGE IF NECESSARY

df["Protein_Class_ID"]=df[COLUMN].map(name_to_id)

df["Original_Class"]=df[COLUMN]

df["PANTHER_Parent_Class"]=df["Protein_Class_ID"].apply(parent_name)

df["PANTHER_Top_Level_Class"]=df["Protein_Class_ID"].apply(top_level)

df["Hierarchy"]=df["Protein_Class_ID"].apply(hierarchy_string)

df["Final_Broad_Class"]=df["Hierarchy"].apply(broad_class)

# ============================================================
# SAVE
# ============================================================

df.to_csv(OUTPUT_FILE,index=False)

print("Done!")

print(df["Final_Broad_Class"].value_counts())

print(f"\nSaved to\n{OUTPUT_FILE}")