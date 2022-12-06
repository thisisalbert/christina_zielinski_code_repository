#Title: .h5ad file of Th17 cells, pre-processing, UMAP, leiden clustering, IL1A expression, differential gene expression
#Figures: Fig. 1a, b, supplementary Fig. 1
#Author: Laurens Lehner

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import seaborn as sns
import collections

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white', figsize=(10, 8), format='png')

path = '/home/laurens/Schreibtisch/Lab/Th17/filtered_feature_bc_matrix'
adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)

#preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.scale(adata, max_value=10)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.4) #0.4 catches most IL1A expressing cells in one cluster

#subsets by cluster
cluster1 = adata[adata.obs['leiden'].isin(['1']),:]
rest = adata[~adata.obs['leiden'].isin(['1']),:]

#add annotation column
cluster1.obs['cluster_identity'] = 'cluster 1'
rest.obs['cluster_identity'] = 'rest'

#merge again
all_cells = cluster1.concatenate(rest)

#subset by cells producing IL1A or not
IL1A = all_cells[all_cells[:,'IL1A'].X>0, :]
noIL1A = all_cells[all_cells[:,'IL1A'].X<=0, :]

#add annotation column
IL1A.obs['contains_IL1A'] = 'IL1A+'
noIL1A.obs['contains_IL1A'] = "IL1A-"

#check how many IL1A+ cells are in cluster 1 and rest
print(collections.Counter(IL1A.obs['cluster_identity']))

#merge again
all_cells = IL1A.concatenate(noIL1A)
all_cells.write_h5ad('all_Th17_cells_no_salt.h5ad')

#all_cells now contains two additional features so we can split by cluster and IL1A expressing cells
#visualize
sc.pl.umap(all_cells, color=['leiden','IL1A','GSDME'], s=200)
sc.pl.umap(all_cells, color=['leiden','IL1A','contains_IL1A','cluster_identity'])

#annotate cells expressing GSDME and cells expressing both GSDME and IL1A
GSDME = all_cells[all_cells[:,'GSDME'].X>0, :]
GSDME_and_IL1A = GSDME[GSDME[:,'IL1A'].X>0, :]
GSDME_and_no_IL1A = GSDME[GSDME[:,'IL1A'].X<0, :]

GSDME_and_IL1A.obs['GSDME_and_IL1A'] = 'GSDME+/IL1A+'
GSDME_and_no_IL1A.obs['GSDME_and_IL1A'] = 'GSDME+/IL1A-'

GSDME =  GSDME_and_IL1A.concatenate(GSDME_and_no_IL1A)

noGSDME = all_cells[all_cells[:,'GSDME'].X<=0, :]
noGSDME_and_IL1A = noGSDME[noGSDME[:,'IL1A'].X>0, :]
noGSDME_and_noIL1A = noGSDME[noGSDME[:,'IL1A'].X<0, :]
noGSDME_and_IL1A.obs['GSDME_and_IL1A'] = 'GSDME-/IL1A+'
noGSDME_and_noIL1A.obs['GSDME_and_IL1A'] = 'GSDME-/IL1A-'
noGSDME = noGSDME_and_IL1A.concatenate(noGSDME_and_noIL1A)
GSDME.obs['contains_GSDME'] = 'GSDME+'
noGSDME.obs['contains_GSDME'] = "GSDME-"

#check how many GSDME+ cells are in cluster 1 and rest
print(collections.Counter(GSDME.obs['cluster_identity']))
#GSDME+ cells expressing IL1A
print(collections.Counter(GSDME.obs['contains_IL1A']))
#GSDME+/IL1A+ cells by cluster
print(collections.Counter(GSDME_and_IL1A.obs['cluster_identity']))

#merge again
all_cells = GSDME.concatenate(noGSDME)
#IL1A+ cells expressing GSDME
IL1A = all_cells[all_cells[:,'IL1A'].X>0, :]
print(collections.Counter(IL1A.obs['contains_GSDME']))

#UMAP with leiden clusters and IL1A expression
sc.pl.umap(all_cells, color=['leiden','IL1A'], s=100,cmap=sns.blend_palette(["lightgrey", sns.xkcd_rgb["black"]], as_cmap=True))

#differentially expressed genes by cluster
sc.tl.rank_genes_groups(all_cells, groupby="leiden")
sc.pl.rank_genes_groups(all_cells, groups="leiden", n_genes=20)