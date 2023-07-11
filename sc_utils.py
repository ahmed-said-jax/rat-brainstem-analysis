from anndata import AnnData
import scanpy as sc
import numpy as np

def cluster_anndata(adata: AnnData,
                    not_highly_variable: list=[],
                    use_existing_neighbors: bool=False,
                    n_neighbors: int=15,
                    n_pcs: int=None,
                    leiden_resolution: float=1,
                    obs_keys_legend_on_data: list=['leiden'],
                    obs_keys_legend_off_data: list=None) -> AnnData:
    '''Perform PCA, compute neighborhood graph, embed with UMAP, and peform Leiden clustering analysis.
    Expects filtered anndata object, not yet normalized or logarithmized.'''

    # Save the raw data in 'raw' layer, normalize, logarithmize, save each in a layer
    adata.layers['raw'] = adata.X.copy()


    sc.pp.normalize_total(adata)
    adata.layers['normalized'] = adata.X.copy()

    sc.pp.log1p(adata)
    adata.layers['log_normed'] = adata.X.copy()

    # Calculate highly variable genes
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger')

    # Exclude certain genes from from PCA, so mark them as not highly variable
    for gene_type in not_highly_variable:
        adata.var.loc[adata.var[gene_type], 'highly_variable'] = False

    # Perform PCA, compute neighborhood graph, embed with UMAP, and cluster with leiden algorithm
    if not use_existing_neighbors:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution)

    if obs_keys_legend_off_data:
        sc.pl.umap(adata, color=obs_keys_legend_off_data)
    
    # Plot results of leiden clustering and the effect of the cells coming from different portions of the brain
    sc.pl.umap(adata, color=obs_keys_legend_on_data, legend_loc='on data')

    return adata

def find_marker_genes(adata: AnnData,
                      filt: bool=False,
                      min_in_group_fraction: float=0,
                      max_out_group_fraction: float=np.inf,
                      n_genes: int=5) -> AnnData:
    
    '''After clustering, find marker genes'''

    # Rank genes according to variability for each cluster
    sc.tl.rank_genes_groups(adata, groupby='leiden')
    ranking_keys = ['rank_genes_groups']
    
    if filt:
        sc.tl.filter_rank_genes_groups(adata, groupby='leiden',
                                       min_in_group_fraction=min_in_group_fraction,
                                       max_out_group_fraction=max_out_group_fraction)
        
        ranking_keys.append('rank_genes_groups_filtered')
    
    # Define URL parameters
    base_url = r'https://panglaodb.se/search.html?query='
    params=r'&species=2&tumor=0&nonadult=0'

    # Some genes are not in panglaodb.se. They were searched, and their correct name was found
    unavailable_genes = {'ENSRNOG00000006030': 'Ptprz1',
                         'H2az1': 'H2afz',
                         'ENSRNOG00000047931': 'Tmsb4x',
                         'Cdc45-1': 'Cdc45',
                         'ENSRNOG00000026336': 'Gas6',
                         'Hmgb2l1-1': 'Hmgb2',
                         'ENSRNOG00000000940': 'Flt1',
                         'ENSRNOG00000019926': 'Ramp1',
                         'Cfap20dc': 'Cfap20'}
    
    # Get the clusters and the rank genes groups dataframe(s)
    clusters = np.unique(adata.obs_vector('leiden'))

    for key in ranking_keys:
        df = sc.get.rank_genes_groups_df(adata, group=clusters, key=key).replace(to_replace=unavailable_genes)

        for cluster in clusters:
            top_genes = df[df['group'] == cluster].nlargest(n=n_genes, columns='scores')['names'].values

            df.loc[df['group'] == cluster, 'panglao_link'] = (
                base_url + '%20AND%20'.join(top_genes) + params
            )

        df.to_pickle(f'{key}.pickle')

    return adata