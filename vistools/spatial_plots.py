import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar

def get_pixel_size_visium(adata, library_id, spot_diameter_real = 55, img='lowres'):
    '''
    Utility function to get size of pixels from visium image in AnnData.
    Useful to plot scalebars with matplotlib_scalebar.ScaleBar

    Params:
    -------
    - adata: AnnData object storing image info in adata.uns[spatial]
    - library_id: string storing ID for image, must be a key in adata.uns[spatial]
    - spot_diameter_real: numeric storing real diameter of spot, in whatever unit you need
        (default is 55, for 55 microns of Visium 10X spots)
    - img: which image to use (default: lowres)
    '''

    ## get scale factor converting original pixel positions (adata.obsm['spatial']) to
    # pixel positions in image
    scalef = adata.uns['spatial'][library_id]['scalefactors']['tissue_{i}_scalef'.format(i=img)]
    ## get spot diameter in image pixels
    spot_diameter_img = adata.uns['spatial'][library_id]['scalefactors']['spot_diameter_fullres'] * scalef
    ## Calculate pixel size
    pixel_size_real = spot_diameter_real/spot_diameter_img
    return(pixel_size_real)

def plot_spatial_scalebar(adata, markers2plot):
    lib_id = [x for x in adata.uns['spatial'].keys()][0]
    pix_size = get_pixel_size_visium(adata=adata, library_id=lib_id)

    fig, ax = plt.subplots()
    sc.pl.spatial(adata, library_id=lib_id, color=markers2plot, use_raw=True,
    size=1.3, img_key='hires',show=False, ax=ax, vmin=0,
    vmax='p99.2', alpha=0.7); # limit color scale at 99.2% quantile of cell abundance

    scalebar = ScaleBar(pix_size, units="um", length_fraction=0.25, frameon=False, location='lower right')
    ax.add_artist(scalebar)
