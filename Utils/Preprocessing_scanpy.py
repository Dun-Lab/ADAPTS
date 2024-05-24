'''
This script includes some helper function for data loading, 
quality control and data normalization on the single cell RNA-seq data 
using scanpy.
'''

import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np

def load_adata(scrna_path, batch_id):
    adata = sc.read_mtx(scrna_path + '/' + batch_id + '_matrix.mtx').T
    adata.var_names = pd.read_csv(scrna_path + '/' + batch_id + '_genes.tsv', header=None, sep='\t')[1]
    adata.obs_names = pd.read_csv(scrna_path + '/' + batch_id + '_barcodes.tsv', header=None, sep='\t')[0]
    adata.obs['batch'] = batch_id
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def basic_qc(adata, quietly=False):

    # Remove duplicated gene_id and cell_id
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Highly expressed genes
    sc.pl.highest_expr_genes(adata, n_top=20, showfliers=False)

    # Annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # QC metrics
    if quietly == False:
        print("Before cell filtering:")
        sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        )

    # Filter outlier cells according to the distribution.
    x = adata.obs['n_genes_by_counts'].values
    x_min = np.max([0, np.median(x)-3*np.std(x)])
    x_max = np.median(x)+3*np.std(x)
    adata = adata[adata.obs.n_genes_by_counts > x_min, :] # adopt to violin plot
    adata = adata[adata.obs.n_genes_by_counts > x_max, :] # adopt to violin plot
    adata = adata[adata.obs.pct_counts_mt < 20, :] # adopt to violin plot
    if quietly == False:
        print("After cell filtering:")
        sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        )
    
    return adata

def Normalized_data(adata):
    # Save the linear scale data
    adata.layers['raw'] = adata.X.copy()

    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e5)
    adata.layers['norm'] = adata.X.copy()

    # Logarithmize the data
    sc.pp.log1p(adata)
    adata.layers['log_norm'] = adata.X.copy()

    # Regress and scale the data
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    adata.layers['scale_data'] = adata.X.copy()
    return adata

def Batch_correct(adata, method='harmony', ):

    # Perform batch correction
    if method == 'raw':
        # Compute distances in the PCA space, and find cell neighbors
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)

    elif method == 'harmony':
        # Batch correction by harmony
        sce.pp.harmony_integrate(adata, 'batch')

        # Compute distances in the PCA space, and find cell neighbors
        sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=10, n_pcs=50)

    elif method == 'scanorama':
        # Batch correction by scanorama
        sce.pp.scanorama_integrate(adata, 'batch', verbose=1)

        # Compute distances in the PCA space, and find cell neighbors
        sc.pp.neighbors(adata, use_rep='X_scanorama', n_neighbors=10, n_pcs=50)
        
    elif method == 'bbknn':
        # Batch correction by BBKNN
        sce.pp.bbknn(adata, 'batch')
    
    elif method == 'mnn':
        # Batch correction by MNN
        hvg = adata.var_names[adata.var['highly_variable']].to_list()
        adata = sce.pp.mnn_correct(adata, var_subset=hvg, batch_key='batch')[0][0]

        # Compute distances in the PCA space, and find cell neighbors
        sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=10, n_pcs=50)

    return adata
    
