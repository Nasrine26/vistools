import scanpy as sc
import numpy as np
import anndata
import pandas as pd
from pathlib import Path
import scipy
from typing import Union, Dict
from matplotlib.image import imread
import json

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

def adata_2_visium(adata: anndata.AnnData,
                   path: Union[str, Path],
                   library_id: str = None,
                  ) -> anndata.AnnData:
    
    r"""
    Creates a visium anndata from an adata.Anndata h5ad annotated count matrix file and 10X-Genomics-formatted visium data directory. 
    The function requires a `spatial` folder where images, scale factors and coordinates can be found. 
    This folder structure is based on the Space Ranger output. 
    Structure of folder should resemble the Space Ranger output for Visium 10X Genomics data.
    
    :param adata: Anndata. The .h5ad annotated data matrix, where spots are identified by their barcodes (rows) and genes by their gene name (columns)
    :param path: Path to directory storing visium files. 
    :param library_id: Identifier for the visium library. 
    
    Inspired from https://github.com/scverse/scanpy/blob/ed3b277b2f498e3cab04c9416aaddf97eec8c3e2/scanpy/readwrite.py#L389
    
    """
    # change path to Path if str
    path = Path(path)
    
    # files needed for spatial
    files = dict(
        tissue_positions_file = path / "spatial/tissue_positions_list.csv",
        scalefactors_json_file = path / "spatial/scalefactors_json.json",
        hires_image = path / "spatial/tissue_hires_image.png",
        lowres_image = path / "spatial/tissue_lowres_image.png"
    )
    
    # create adata.uns['spatial'][<library_id>] slot, where spatial info referring to tissue image is stored
    # follow scanpy structure 
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][library_id] = dict()
    
    # create key slot for images
    adata.uns["spatial"][library_id]["images"] = dict()
    
    # read in images
    for res in ["hires", "lowres"]:
        try:
            img = imread(str(files[f"{res}_image"]))
            adata.uns["spatial"][library_id]["images"][res] = img
            
        except Exception:
            raise OSError(f"Could not find '{res}_image'")
            
    # store image path in meta aka path to the high-resolution tissue image
    adata.uns["spatial"][library_id]["metadata"] = dict()
    source_image_path = str(files["hires_image"])
    adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = source_image_path
    
    
    # read scale factors from json file
    scalefactors = json.loads(files["scalefactors_json_file"].read_bytes())
    adata.uns["spatial"][library_id]["scalefactors"] = scalefactors
    
    # read coordinates from csv file
    positions = pd.read_csv(files['tissue_positions_file'], header=None)
    
    
    positions.columns = [
        'barcode',
        'in_tissue',
        'array_row',
        'array_col',
        'pxl_col_in_fullres',
        'pxl_row_in_fullres'
    ]
    
    positions.index = positions['barcode']

    adata.obs = adata.obs.merge(positions, how="left", left_index=True, right_index=True)

    adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()

    adata.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)
    
    return adata



