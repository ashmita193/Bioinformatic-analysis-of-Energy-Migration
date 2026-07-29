# Bioinformatic-analysis-of-Energy-Migration
Data Analysis workflow of Diffusion length data. 

Scripts

This repository contains Python scripts for preprocessing, statistical analysis, and visualization of protein excitation energy migration (diffusion length) data. Each script performs a specific stage of the computational workflow.

1. residue_count.py

Retrieves protein sequences from UniProt and calculates the aromatic amino acid composition for each protein.

Functions
Downloads protein sequences using UniProt accession IDs
Counts Tryptophan (W), Tyrosine (Y), and Phenylalanine (F) residues
Calculates total aromatic residue content
Generates a sequence composition dataset for downstream analyses

Input

CSV file containing UniProt accession IDs

Output

CSV containing protein length and aromatic residue counts

2. ekthapanther.py

Annotates proteins using the PANTHER Protein Classification database and assigns functional protein classes.

Functions
Maps UniProt IDs to PANTHER protein classes
Retrieves functional annotations
Standardizes class names
Generates a protein classification dataset

Input

Protein accession IDs

Output

Protein functional classification table

3. density_analysis.py

Performs exploratory analysis of diffusion length distributions.

Functions
Computes kernel density estimates (KDE)
Generates diffusion length distribution plots
Examines distribution shape and variability

Input

Diffusion length dataset

Output

Density plots
Distribution statistics

4. for_science.py

Performs statistical comparisons between protein classes.

Functions
Organizes proteins by functional class
Performs Kruskal–Wallis tests
Conducts Dunn's post hoc multiple comparisons
Computes adjusted p-values
Generates statistical summary tables

Input

Protein classification dataset
Diffusion length dataset

Output

Statistical test results
Pairwise comparison tables

5. PCA.py

Performs Principal Component Analysis (PCA) to identify major sources of variation within the dataset.

Functions
Standardizes numerical variables
Computes principal components
Calculates explained variance
Generates PCA score plots
Produces loading plots for feature interpretation

Input

Combined protein feature dataset

Output

PCA plots
Principal component loadings
Explained variance statistics

6. broad_class_violin.py

Visualizes diffusion length distributions across broad protein functional classes for FRET data.

Functions
Generates violin plots
Displays medians and distribution density
Compares diffusion lengths among protein classes

Input

Protein classification dataset
Diffusion length dataset

Output

Violin plots

7. broad_class_violin_det.py

Visualizes diffusion length distributions across broad protein functional classes for DET data.

Functions
Applies custom ordering of protein classes
Uses publication-ready colour schemes
Generates high-resolution figures
Includes statistical annotations where applicable

Input

Classified diffusion length dataset

Output

Publication-quality violin plots

8. heatmap.py

Generates heatmaps for correlation and comparative analyses.

Functions
Computes correlation matrices
Visualizes relationships among variables
Produces publication-quality heatmaps

Input

Numerical feature dataset

Output

Correlation heatmaps
Comparative heatmaps


Recommended Execution Order
1. residue_count.py
        │
        ▼
2. ekthapanther.py
        │
        ▼
3. Merge datasets
        │
        ▼
4. density_analysis.py
        │
        ▼
5. for_science.py
        │
        ▼
6. PCA.py
        │
        ▼
7. broad_class_violin.py
        │
        ▼
8. broad_class_violin_det.py
        │
        ▼
9. heatmap.py

Dependencies

The scripts primarily use the following Python libraries:

pandas
numpy
matplotlib
scipy
scikit-learn
seaborn
requests
BioPython
statsmodels
scikit-posthocs
Outputs

The workflow generates:

Protein aromatic residue composition datasets
Protein functional classification tables
Principal Component Analysis (PCA) visualizations
Kernel density estimation (KDE) plots
Violin plots of diffusion length distributions
Correlation heatmaps
Kruskal–Wallis and Dunn's post hoc statistical results
Publication-ready figures suitable for scientific manuscripts
