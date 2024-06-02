'''
This script is used to load the single cell RNA-seq data using scanpy.
'''
import pandas as pd
import scanpy as sc

def load_adata(scrna_path, batch_id):
    adata = sc.read_mtx(scrna_path + '/' + batch_id + '_matrix.mtx').T
    adata.var_names = pd.read_csv(scrna_path + '/' + batch_id + '_genes.tsv', header=None, sep='\t')[1]
    adata.obs_names = pd.read_csv(scrna_path + '/' + batch_id + '_barcodes.tsv', header=None, sep='\t')[0]
    adata.obs['batch'] = batch_id
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

