import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import sklearn.linear_model 
from sklearn import metrics
from sklearn.metrics import mean_squared_error
import anndata
import scanpy as sc
import scipy
from patsy import dmatrix

def contraster(dds_coldata, design_formula, group1, group2, weighted=False):
    """
    dds_coldata: pandas DataFrame (equivalent to colData(dds))
    design_formula: string formula (e.g. '~ condition + batch')
    group1: list of lists, where each list contains a column name followed by values
    group2: same as group1
    weighted: boolean flag for including duplicate rows
    """
    
    # Create model matrix
    mod_mat = dmatrix(design_formula, dds_coldata, return_type='dataframe')
    
    # Logical index for group1
    grp1_rows = [dds_coldata[g[0]].isin(g[1:]) for g in group1]
    grp1_mask = np.logical_and.reduce(grp1_rows)
    
    # Logical index for group2
    grp2_rows = [dds_coldata[g[0]].isin(g[1:]) for g in group2]
    grp2_mask = np.logical_and.reduce(grp2_rows)
    
    # Subset model matrices
    mod_mat1 = mod_mat.loc[grp1_mask]
    mod_mat2 = mod_mat.loc[grp2_mask]
    
    if not weighted:
        mod_mat1 = mod_mat1.drop_duplicates()
        mod_mat2 = mod_mat2.drop_duplicates()
    
    # Return the difference in column means
    return (mod_mat1.mean() - mod_mat2.mean()).to_numpy()

def nd(arr):
    """
    Function to transform numpy matrix to nd array.
    """
    return np.asarray(arr).reshape(-1)
    
def sum_by(adata: anndata.AnnData, col: str) -> anndata.AnnData:
    adata.strings_to_categoricals()
    assert pd.api.types.is_categorical_dtype(adata.obs[col])

    cat = adata.obs[col].values
    indicator = scipy.sparse.coo_matrix(
        (
            np.broadcast_to(True, adata.n_obs),
            (cat.codes, np.arange(adata.n_obs))
        ),
        shape=(len(cat.categories), adata.n_obs),
    )

    return anndata.AnnData(
        indicator @ adata.X,
        var=adata.var,
        obs=pd.DataFrame(index=cat.categories)
    )
    
def genes_in_original(gene_list, og_supp_table =2, t2go_path = 't2go.csv', og_supp_table_path = 'metadata_csvs/',tx = False):
    t2go = pd.read_csv(t2go_path)
    if og_supp_table ==2:
        og_de = pd.read_csv(og_supp_table_path+f'og_supp{og_supp_table}.csv', skiprows = 18)
        og_de = og_de[og_de.assay == 'TRNSCRPT']
        table_lookup = pd.merge(og_de, t2go, left_on = 'feature_ID', right_on = 'ensembl_gene_id', how = 'left')
    elif og_supp_table == 3:
        og_de = pd.read_csv(og_supp_table_path+f'og_supp{og_supp_table}.csv', skipfooter=8)
        table_lookup = pd.merge(og_de, t2go, left_on = 'gene_symbol', right_on = 'X3', how = 'left')
    else:
        og_de = pd.read_csv(og_supp_table_path+f'og_supp{og_supp_table}.csv', skiprows = 10)
        og_de = og_de[og_de.assays == 'TRNSCRPT']
        table_lookup = pd.merge(og_de, t2go, left_on = 'gene_symbol', right_on = 'X3', how = 'left')

    if tx:
        return([x in table_lookup.ensembl_transcript_id.unique().tolist() for x in gene_list])
    else:
        return([x in table_lookup.ensembl_gene_id.unique().tolist() for x in gene_list])

def get_CCC(real, pred):
    cor = np.corrcoef(real, pred)[0][1]
    # Means
    mean_true = np.mean(real)
    mean_pred = np.mean(pred)
    # Population variances
    var_true = np.var(real)
    var_pred = np.var(pred)
    # Population standard deviations
    sd_true = np.std(real)
    sd_pred = np.std(pred)
    # Calculate CCC
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred)**2
    ccc = numerator / denominator
    return ccc

def run_and_eval(X, y, z, sex = 1, reg = 'ridge', state = 0):
    if sex == 1:
        train_filter = y[['pid', 'time']].groupby('pid').agg(lambda x: x.value_counts().index[0])
        train_pids = train_filter.groupby(['time']).apply(lambda x:x.sample(1+np.round(len(x)/3).astype('int'))).index.get_level_values(1).tolist()
    else:
        train_filter = y[['pid', 'time', 'sex']].groupby('pid').agg(lambda x: x.value_counts().index[0])
        train_pids = train_filter.groupby(['time', 'sex']).apply(lambda x:x.sample(1+np.round(len(x)/3).astype('int'))).index.get_level_values(2).tolist()
    ref_mask = [x in train_pids for x in y.pid]
    X_train = X[ref_mask]
    X_test = X[list(~np.array(ref_mask))]
    y_train = y[ref_mask]
    y_test = y[list(~np.array(ref_mask))]
    if reg == 'ridge':
        model = sklearn.linear_model.Ridge(random_state=0)
    else:
        model = sklearn.linear_model.Lasso(random_state=0)
    model.fit(X_train, y_train.time)
    y_pred = model.predict(X_test)
    return get_CCC(y_test.time.astype('int'), y_pred), ref_mask

def load_annotated_omic(adata_path, omic):
    adata = anndata.read_h5ad(adata_path)
    if omic in ['RNA']:
        adata.layers["counts"] = adata.X.copy()
        adata.obs.time = adata.obs.time.astype('int')
        adata.var = pd.DataFrame({'X': adata.var.gene_id})
        adata.var['omic'] = 'RNA'
    else:
        adata.X = np.nan_to_num(adata.X, nan=0)
        adata.var['omic'] = omic
        #adata.var.index = adata.var.X
        adata.obs.index = adata.obs.pid.astype('int').astype('str')
    adata.layers['counts'] = adata.X
    return adata

def data_reshaper_omic(a_data, tiss_list, ngenes = 90, mode = 'all',omic='RNA'):
    adata = a_data.copy()        
    start = True
    if mode == 'hvg':
        sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=ngenes,
                    subset=True,
                    layer="counts",
                    flavor="seurat_v3",
                )
    for tiss in tiss_list:
        # tissues missing metadata removed
        if (tiss in ['OVARY', 'TESTES', 'VENACV']):
            continue
        print(tiss)
        if start:
            big_data = adata[adata.obs.tissue==tiss].copy()
            big_data.var_names = big_data.var_names + f'_{tiss}'
            big_data.var['tissue'] = tiss
            big_data.obs_names = big_data.obs.nid.astype('str')
            start=False
        else:
            curr_data = adata[adata.obs.tissue==tiss].copy()
            curr_data.var_names = curr_data.var_names + f'_{tiss}'
            curr_data.var['tissue'] = tiss
            curr_data.obs_names = curr_data.obs.nid.astype('str')
            big_data = anndata.concat([big_data,curr_data],axis=1, merge='same')
    big_data = big_data[:,big_data.layers['counts'].sum(axis=0)>0]
    slim_pheno = a_data.obs.drop('tissue',axis=1).reset_index(drop=True)
    if omic in ['ATAC', 'METHYL']:
        big_data.obs = pd.merge(big_data.obs.reset_index(drop=True), slim_pheno.drop_duplicates('nid'), on=['nid'], how='left', suffixes=('', '_DROP')).filter(regex='^(?!.*_DROP)').set_index('nid')
    else:
        big_data.obs = pd.merge(big_data.obs.reset_index(names='nid'), slim_pheno.drop_duplicates('nid'), on=['nid'], how='left', suffixes=('', '_DROP')).filter(regex='^(?!.*_DROP)').set_index('nid')
    # filter for na
    big_data = big_data[~big_data.obs.time.isna()]
    z = big_data.var
    y = big_data.obs
    y.reset_index(inplace=True)
    X = pd.DataFrame(big_data.layers['counts'],
                     columns = big_data.var_names)
    return(X, y, z)

def lts(pvalue, weight):
    return scipy.stats.gamma.ppf(1-pvalue, a=weight/2, scale=2)
    
def lancaster(pvalues, weights):
    if len(weights) != len(pvalues):
        raise ValueError("Length of weights not equal to length of pvalues")
    
    # Remove NaN values
    valid_mask = ~np.isnan(pvalues)
    pvalues = np.array(pvalues)[valid_mask]
    weights = np.array(weights)[valid_mask]
    
    # Apply weight condition
    valid_weight_mask = (weights > 0) | np.isnan(weights)
    pvalues = pvalues[valid_weight_mask]
    weights = weights[valid_weight_mask]
    
    if len(pvalues) == 0:
        return np.nan
    if len(pvalues) == 1:
        return pvalues[0]
    
    # Check for extreme p-values
    if np.any(pvalues < 9.99988867182683e-320):
        print("Warning: Extreme p-values around and below 10e-320 will produce a p-value of 0. Replace extreme p-values with 10e-320 to obtain an upper bound for the aggregated p-value.")
    
    # Calculate aggregated t-values
    t = np.array([lts(p, w) for p, w in zip(pvalues, weights)])
    
    # Sum the t-values and calculate the chi-squared p-value
    t_sum = np.sum(t)
    p = scipy.stats.chi2.sf(t_sum, np.sum(weights))  # sf is the survival function (1 - cdf)
    
    return p