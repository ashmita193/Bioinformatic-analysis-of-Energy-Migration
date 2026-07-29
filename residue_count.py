import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ==========================================================
# INPUT / OUTPUT
# ==========================================================
INPUT_CSV = r"C:\Users\91892\kalralabintern\pr_fams\proteinids.csv"
OUTPUT_CSV = r"C:\Users\91892\kalralabintern\protein_with_aromatic_residue_counts.csv"

# Number of threads
N_THREADS = 20

# ==========================================================
# READ INPUT
# ==========================================================
df = pd.read_csv(INPUT_CSV)

# ==========================================================
# CREATE A SESSION WITH RETRIES
# ==========================================================
session = requests.Session()

retry = Retry(
    total=5,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
)

adapter = HTTPAdapter(max_retries=retry)
session.mount("https://", adapter)

# ==========================================================
# FUNCTION TO DOWNLOAD SEQUENCE AND COUNT RESIDUES
# ==========================================================
def count_residues(uniprot_id):

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    try:
        r = session.get(url, timeout=20)

        if r.status_code != 200:
            return {
                "Uniprot_ID": uniprot_id,
                "Length": None,
                "Tryptophan(W)": None,
                "Tyrosine(Y)": None,
                "Phenylalanine(F)": None,
            }

        fasta = r.text

        sequence = "".join(
            line.strip()
            for line in fasta.splitlines()
            if not line.startswith(">")
        )

        return {
            "Uniprot_ID": uniprot_id,
            "Length": len(sequence),
            "Tryptophan(W)": sequence.count("W"),
            "Tyrosine(Y)": sequence.count("Y"),
            "Phenylalanine(F)": sequence.count("F"),
        }

    except Exception:
        return {
            "Uniprot_ID": uniprot_id,
            "Length": None,
            "Tryptophan(W)": None,
            "Tyrosine(Y)": None,
            "Phenylalanine(F)": None,
        }

# ==========================================================
# MULTITHREADED DOWNLOAD
# ==========================================================
results = []

with ThreadPoolExecutor(max_workers=N_THREADS) as executor:

    futures = {
        executor.submit(count_residues, uid): uid
        for uid in df["Uniprot_ID"]
    }

    for i, future in enumerate(as_completed(futures), 1):
        results.append(future.result())

        if i % 100 == 0:
            print(f"{i} proteins processed")

# ==========================================================
# MERGE RESULTS
# ==========================================================
result_df = pd.DataFrame(results)

final_df = df.merge(result_df, on="Uniprot_ID", how="left")

# ==========================================================
# SAVE
# ==========================================================
final_df.to_csv(OUTPUT_CSV, index=False)

print(f"\nSaved to {OUTPUT_CSV}")