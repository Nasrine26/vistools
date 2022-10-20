import scanpy as sc
import numpy as np

def read_and_qc(path, count_file_prefix='', sample_name=None):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics

    :param sample_name: Name of the sample
    :param path: path to data
    :param count_file_prefix: prefix in front of count file name filtered_feature_bc_matrix.h5
    """

    adata = sc.read_visium(path,
                           count_file=f'{count_file_prefix}filtered_feature_bc_matrix.h5', load_images=True)

    print('Sample ', (list(adata.uns['spatial'].keys())[0]))
    adata.obs['Sample'] = list(adata.uns['spatial'].keys())[0]
    #adata.obs['sample'] = sample_name

    #adata.var['SYMBOL'] = adata.var_names
    #adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    #adata.var_names = adata.var['ENSEMBL']
    #adata.var.drop(columns='ENSEMBL', inplace=True)

    # since we need unique gene names (or it would actually be better to have ensembl ids), we can make them unique
    # Otherwise, error occurs (intersect, find shared genes and subset both anndata and reference signatures)
    adata.var_names_make_unique()

    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()

    # find mitochondria-encoded (MT) genes
    adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var_names]
    adata.var['ribo'] = [gene.startswith(('RPS','RPL')) for gene in adata.var_names]
    #adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=False, qc_vars=['mt','ribo'])
    adata.X = csr_matrix(adata.X)

    # add sample name to obs names
    adata.obs["Sample"] = [str(i) for i in adata.obs['Sample']]
    adata.obs_names = adata.obs["Sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata

def select_slide(adata, s, batch_key="sample"):
    r"""This function selects the data for one slide from the spatial anndata object (multiple samples).

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param batch_key: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[batch_key].isin([s]), :].copy()
    s_keys = list(slide.uns["spatial"].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]

    slide.uns["spatial"] = {s_spatial: slide.uns["spatial"][s_spatial]}

    return slide
