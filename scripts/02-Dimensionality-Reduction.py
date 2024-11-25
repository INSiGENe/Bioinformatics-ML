import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Load data
adata = sc.datasets.pbmc3k()

# Calculate mitochondrial genes percentage
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Quality control
# Filter cells based on number of genes
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# Filter cells based on mitochondrial percentage
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Feature selection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# Visualize QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Plot relationships between QC metrics
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

# Compare different numbers of PCs for variance explained
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(adata.uns['pca']['variance_ratio'])
plt.xlabel('PC')
plt.ylabel('Variance ratio')
plt.title('Explained variance ratio by PC')

plt.subplot(1, 2, 2)
plt.plot(np.cumsum(adata.uns['pca']['variance_ratio']))
plt.xlabel('Number of PCs')
plt.ylabel('Cumulative variance ratio')
plt.title('Cumulative explained variance ratio')
plt.tight_layout()
plt.show()

# t-SNE with different perplexities
perplexities = [5, 30, 50]
fig, axes = plt.subplots(1, len(perplexities), figsize=(15, 5))

for idx, perp in enumerate(perplexities):
    # Create a copy of the data for each perplexity value
    adata_copy = adata.copy()
    # Run t-SNE
    sc.tl.tsne(adata_copy, use_rep='X_pca', perplexity=perp)
    # Plot
    sc.pl.tsne(adata_copy, color='CST3', 
               title=f'perplexity={perp}', 
               ax=axes[idx], 
               show=False)

plt.tight_layout()
plt.show()

# UMAP with different n_neighbors
n_neighbors_list = [5, 15, 30]
fig, axes = plt.subplots(1, len(n_neighbors_list), figsize=(15, 5))

for idx, n_neigh in enumerate(n_neighbors_list):
    # Create a copy of the data for each n_neighbors value
    adata_copy = adata.copy()
    # Calculate neighbors and run UMAP
    sc.pp.neighbors(adata_copy, n_neighbors=n_neigh, n_pcs=30)
    sc.tl.umap(adata_copy)
    # Plot
    sc.pl.umap(adata_copy, color='CST3',
               title=f'n_neighbors={n_neigh}',
               ax=axes[idx],
               show=False)

plt.tight_layout()
plt.show()

# Final combined visualization with default parameters
# Run t-SNE and UMAP with default parameters
sc.tl.tsne(adata, use_rep='X_pca')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

# Create combined plot
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

sc.pl.pca(adata, color='CST3', ax=ax1, show=False)
ax1.set_title('PCA')

sc.pl.tsne(adata, color='CST3', ax=ax2, show=False)
ax2.set_title('t-SNE')

sc.pl.umap(adata, color='CST3', ax=ax3, show=False)
ax3.set_title('UMAP')

plt.tight_layout()
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print(f"Final number of cells: {adata.n_obs}")
print(f"Final number of genes: {adata.n_vars}")
print(f"Number of PCs used: {adata.uns['pca']['variance_ratio'].shape[0]}")
