#!/usr/bin/python3

from pathlib import Path
import scipy
import numpy as np
import os
import anndata as ad
import pandas as pd

output = 'boot'
files = ['matrix.abundance.mtx', 'matrix.abundance.tpm.mtx', 'matrix.abundance.gene.mtx', 'matrix.abundance.gene.tpm.mtx']
for x in files:
    start = True
    pseudo_bar = list()
    for result in Path(f'./').glob('count_matrices_SRR*'):
        pseudo_bar.append(str(result)[-11:])
        print(result)
        if start:
            big_data = scipy.io.mmread(result/f'quants_unfiltered_{output}'/x)
            start = False
        else:
            big_data = scipy.sparse.vstack([big_data, scipy.io.mmread(result/f'quants_unfiltered_{output}'/x)])
    scipy.io.mmwrite(f'{output}_{x}', big_data)
np.savetxt(f'{output}_pseudo_bar.txt', pseudo_bar, delimiter="\n", fmt="%s")

for x in ['matrix.abundance.tpm.mtx', 'matrix.abundance.gene.tpm.mtx']:
    big_data = scipy.io.mmread(f'{output}_{x}').tocsr()
    transcripts = list(pd.read_csv(result/f'quants_unfiltered_{output}'/'transcripts.txt',sep='\t', header=None)[0])
    genes = list(pd.read_csv(result/f'quants_unfiltered_{output}'/'genes.txt',sep='\t', header=None)[0])
    barcodes = list(pd.read_csv(f'{output}_pseudo_bar.txt',sep='\t', header=None)[0])
    if x == 'matrix.abundance.tpm.mtx':
        adata = ad.AnnData(big_data, obs=pd.DataFrame(barcodes, columns=['barcode']), var=pd.DataFrame(transcripts, columns=['transcript']))
        adata.write(f'{output}_transcripts_tpm.h5ad', compression="gzip")
    else:
        adata = ad.AnnData(big_data, obs=pd.DataFrame(barcodes, columns=['barcode']), var=pd.DataFrame(genes, columns=['genes']))
        adata.write(f'{output}_genes_tpm.h5ad', compression="gzip")

for x in ['matrix.abundance.mtx', 'matrix.abundance.gene.mtx']:
    big_data = scipy.io.mmread(f'{output}_{x}').tocsr()
    transcripts = list(pd.read_csv(result/f'quants_unfiltered_{output}'/'transcripts.txt',sep='\t', header=None)[0])
    genes = list(pd.read_csv(result/f'quants_unfiltered_{output}'/'genes.txt',sep='\t', header=None)[0])
    barcodes = list(pd.read_csv(f'{output}_pseudo_bar.txt',sep='\t', header=None)[0])
    if x == 'matrix.abundance.mtx':
        adata = ad.AnnData(big_data, obs=pd.DataFrame(barcodes, columns=['barcode']), var=pd.DataFrame(transcripts, columns=['transcript']))
        adata.write(f'{output}_transcripts.h5ad', compression="gzip")
    else:
        adata = ad.AnnData(big_data, obs=pd.DataFrame(barcodes, columns=['barcode']), var=pd.DataFrame(genes, columns=['genes']))
        adata.write(f'{output}_genes.h5ad', compression="gzip")
