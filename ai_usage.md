# AI Use Log

Tool/model & version: ChatGPT

Prompt: Load and parse the GEO series matrix file GSE11877_series_matrix.txt
Given code: import pandas as pd
import numpy as np
import re

file = "GSE11877_series_matrix.txt"

# Read entire file as text lines
with open(file, "r") as f:
    lines = f.readlines()

# Identify where the data table begins and ends
start = lines.index("!series_matrix_table_begin\n") + 1
end = lines.index("!series_matrix_table_end\n")

# Read expression matrix between those markers
expr = pd.read_csv(
    file,
    sep="\t",
    skiprows=start,
    nrows=end - start
)

expr.head()
What I changed before committing: N/A

Prompt: Extract sample metadata for tissue type
Given code: # Pull the tissue-type annotation line
tissue_line = [l for l in lines 
               if l.startswith("!Sample_characteristics_ch1\t\"tissue type")][0]

# Split columns (remove the first token: !Sample_characteristics_ch1)
tissue_labels = tissue_line.split("\t")[1:]

# Clean up quotes and "tissue type :" prefix
tissue_labels = [
    x.replace('"', '').replace("tissue type :", "").strip()
    for x in tissue_labels
]

# Simplify to bone marrow vs peripheral blood labels
tissue_simple = [
    "bone_marrow" if "bone marrow" in x.lower() else "peripheral_blood"
    for x in tissue_labels
]

len(tissue_simple)
What I changed before committing: N/A

Prompt: Load and parse the GEO series matrix file GSE11877_series_matrix.txt
Given code: import pandas as pd
import numpy as np
import re

file = "GSE11877_series_matrix.txt"

with open(file, "r") as f:
    lines = f.readlines()

start = lines.index("!series_matrix_table_begin\n") + 1
end = lines.index("!series_matrix_table_end\n")

expr = pd.read_csv(
    file,
    sep="\t",
    skiprows=start,
    nrows=end - start
)

expr.head()
What I changed before committing: N/A

Prompt: Clean and format the expression matrix
Given code: # First column = probe IDs
expr = expr.rename(columns={expr.columns[0]: "ID_REF"})
expr = expr.set_index("ID_REF")

# Convert all remaining entries to numeric
expr = expr.apply(pd.to_numeric, errors="coerce")

# Remove probes with all missing values
expr = expr.dropna(how="all")

# Build phenotype (tissue) DataFrame
pheno = pd.DataFrame({"tissue": tissue_simple})
pheno.index = expr.columns  # ensure sample order matches columns
What I changed before committing: N/A

Prompt: Perform a differential expression analysis between bone marrow and peripheral blood samples
Given code: from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Define sample groups
bm_samples = pheno[pheno["tissue"] == "bone_marrow"].index
pb_samples = pheno[pheno["tissue"] == "peripheral_blood"].index

bm = expr[bm_samples]
pb = expr[pb_samples]

# Compute log2 fold change (bone marrow vs peripheral blood)
logFC = np.log2(bm.mean(axis=1) + 1) - np.log2(pb.mean(axis=1) + 1)

# Two-sample t-test for each probe
pvals = ttest_ind(bm.T, pb.T, nan_policy="omit").pvalue

# FDR correction (Benjamini–Hochberg)
fdr = multipletests(pvals, method="fdr_bh")[1]

# Collect results into a DataFrame
DE = pd.DataFrame({
    "logFC": logFC,
    "pval": pvals,
    "FDR": fdr
})

DE_sorted = DE.sort_values("FDR")
DE_sorted.head()
What I changed before committing: N/A

Prompt: Generate PCA, volcano, and heatmap plots
Given code: from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# PCA
pca = PCA(n_components=2)
pcs = pca.fit_transform(expr.T)

pca_df = pd.DataFrame({
    "PC1": pcs[:, 0],
    "PC2": pcs[:, 1],
    "tissue": tissue_simple
})

plt.figure(figsize=(7, 6))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="tissue")
plt.title("PCA: Bone Marrow vs Peripheral Blood")
plt.show()

# Volcano plot
plt.figure(figsize=(7, 6))
sns.scatterplot(
    x=DE["logFC"],
    y=-np.log10(DE["pval"]),
    hue=DE["FDR"] < 0.05,
    legend=False
)
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 p-value")
plt.title("Volcano Plot")
plt.show()

# Heatmap of top 30 differentially expressed genes
top_genes = DE_sorted.head(30).index

plt.figure(figsize=(10, 8))
sns.heatmap(expr.loc[top_genes], cmap="viridis")
plt.title("Top 30 Differentially Expressed Genes")
plt.show()

# Save DE table
DE_sorted.to_csv("DE_BM_vs_PB_GSE11877.csv")
What I changed before committing: N/A
