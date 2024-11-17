## This script was used to perform different types of feature selection.   

salloc -A def-sfarhan --time=0-8 -c 1 --mem=50g
module load StdEnv/2020 
module load r/4.2.2 
R


salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

###################################
# Adjustable parameters
###################################
# optimal PCs for PCA
n_pc = 50
# optimal ranks for NMF
n_rank = 24
# optimal topics for ETM
n_topics = 50 


par_test_anndata = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/test_anndata_base.h5ad'
par_train_anndata = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/train_anndata_base.h5ad'
par_m1_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m1_HVGS_train.h5ad'
par_m1_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m1_HVGS_test.h5ad'
par_m2_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m2_PCA_train.h5ad' 
par_m2_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m2_PCA_test.h5ad' 
par_m3_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m3_NMF_factors_train.h5ad'
par_m3_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m3_NMF_factors_test.h5ad'
par_m4_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m4_ETM_train.h5ad'
par_m4_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m4_ETM_test.h5ad'

par_top_genes_topic = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/topics_top200genes.csv'
par_top_genes_factor = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/factors_top_genes.csv'


###################################
# load library
###################################
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, vstack
from numpy import inf
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
import anndata as ad
from scETM import scETM, UnsupervisedTrainer, evaluate, prepare_for_transfer
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

###################################
# scanpy settings
###################################
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


###################################
# Load anndata
###################################
train_anndata = sc.read_h5ad(par_train_anndata)
test_anndata =  sc.read_h5ad(par_test_anndata)


###################################
# Define preprocessing functions
###################################
## Annotate MT genes
def annotate_mt_genes(adata):
  adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
  sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

## Normalize data
def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    return adata

## Logarithmize data
def logarithmize(adata):
    sc.pp.log1p(adata)
    return adata

## Identify HGVs
def HVGs(adata):
    '''Find out the HVGs (without specific number of it) and store in the variable_genes parameter'''
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    return adata

## Select HGVs  
def select_HVGs(adata):
    adata = adata[:, adata.var.highly_variable]  
    print("dimension check:" + str(adata))
    return adata


######################################################################
# Preprocessing
# We are going to preprocess the  new_anndata object
# There will be 4 preprocessing methods
#   1. HVGs 
#   2. PCA
#   3. NMF
#   4. Topics
######################################################################

##################################################################################################################################################################################################################
# Method 1:
##################################################################################################################################################################################################################
## copy the data object
m1_5000_genes_train = train_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m1_5000_genes_train)

## Normalize data
m1_5000_genes_train = normalize(m1_5000_genes_train)

## Logarithmize data
m1_5000_genes_train = logarithmize(m1_5000_genes_train)

## Identify highly variable genes
m1_5000_genes_train = HVGs(m1_5000_genes_train)

## Getting the 5000 HVGs and store in dataframe
print('size before HVGs: ' + str(m1_5000_genes_train.X.shape) + str(m1_5000_genes_train.obs.shape))
m1_5000_genes_train = select_HVGs(m1_5000_genes_train)
print('size after HVGs: ' + str(m1_5000_genes_train.X.shape) + str(m1_5000_genes_train.obs.shape))

## Rename disease ontology
m1_5000_genes_train.obs['disease__ontology_label'] = m1_5000_genes_train.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## Create new anndata object
a_train = sc.AnnData(X = m1_5000_genes_train.X.copy(),
  obs = m1_5000_genes_train.obs.copy(),
  var = m1_5000_genes_train.var.copy(),
  uns = m1_5000_genes_train.uns.copy(),
  obsm = m1_5000_genes_train.obsm.copy(),
  varm = m1_5000_genes_train.varm.copy(),
  layers = m1_5000_genes_train.layers.copy(),
  raw = m1_5000_genes_train.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = m1_5000_genes_train.obsp.copy(),
  varp = m1_5000_genes_train.varp
  )

a_train.__dict__['_raw'].__dict__['_var'] = m1_5000_genes_train.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 1 data object
a_train.write(par_m1_train)


## pull out the same highly variable genes on test dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## copy the data object
m1_5000_genes_test = test_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m1_5000_genes_test)

## Normalize data
m1_5000_genes_test = normalize(m1_5000_genes_test)

## Logarithmize data
m1_5000_genes_test = logarithmize(m1_5000_genes_test)

## only keep HGV gnes from train datasets
m1_5000_genes_test = m1_5000_genes_test[:,m1_5000_genes_train.var_names]

## Rename disease ontology
m1_5000_genes_test.obs['disease__ontology_label'] = m1_5000_genes_test.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## Create new anndata object
a_test = sc.AnnData(X = m1_5000_genes_test.X.copy(),
  obs = m1_5000_genes_test.obs.copy(),
  var = m1_5000_genes_test.var.copy(),
  uns = m1_5000_genes_test.uns.copy(),
  obsm = m1_5000_genes_test.obsm.copy(),
  varm = m1_5000_genes_test.varm.copy(),
  layers = m1_5000_genes_test.layers.copy(),
  raw = m1_5000_genes_test.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = m1_5000_genes_test.obsp.copy(),
  varp = m1_5000_genes_test.varp
  )

a_test.__dict__['_raw'].__dict__['_var'] = m1_5000_genes_test.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 1 data object
a_test.write(par_m1_test)


############################################################################################################################################
# Method 2:
############################################################################################################################################
## PCA model on train data
m2_50_PCA_train = train_anndata.copy()

## Annotate data object
annotate_mt_genes(m2_50_PCA_train)

## Normalize the data
m2_50_PCA_train = normalize(m2_50_PCA_train) 

## Logarithmize the data
m2_50_PCA_train = logarithmize(m2_50_PCA_train)

## PCA model
model = PCA(n_components=n_pc, random_state=0)

## fit the model
model.fit(m2_50_PCA_train.X)
P = model.transform(m2_50_PCA_train.X)

type(P)

## Rename disease ontology
m2_50_PCA_train.obs['disease__ontology_label'] = m2_50_PCA_train.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## check
new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(P.shape[1])])

## Create new anndata object for method 2
b_train = sc.AnnData(X = P.copy(),
  obs = m2_50_PCA_train.obs.copy(),
  var = new_var,
  uns = m2_50_PCA_train.uns.copy(),
  obsm = m2_50_PCA_train.obsm.copy(),
  varm = m2_50_PCA_train.varm.copy(),
  layers = m2_50_PCA_train.layers.copy(),
  raw = m2_50_PCA_train.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = m2_50_PCA_train.obsp.copy(),
  varp = m2_50_PCA_train.varp
  )
b_train.__dict__['_raw'].__dict__['_var'] = m2_50_PCA_train.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 2 data object
b_train.write(par_m2_train)

## Apply PCA model to test set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PCA model on train data
m2_50_PCA_test = test_anndata.copy()

## Annotate data object
annotate_mt_genes(m2_50_PCA_test)

## Normalize the data
m2_50_PCA_test = normalize(m2_50_PCA_test) 

## Logarithmize the data
m2_50_PCA_test = logarithmize(m2_50_PCA_test)

P = model.transform(m2_50_PCA_test.X)

type(P)

## Rename disease ontology
m2_50_PCA_test.obs['disease__ontology_label'] = m2_50_PCA_test.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## check
new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(P.shape[1])])


## Create new anndata object for method 2
b_test = sc.AnnData(X = P.copy(),
  obs = m2_50_PCA_test.obs.copy(),
  var = new_var,
  uns = m2_50_PCA_test.uns.copy(),
  obsm = m2_50_PCA_test.obsm.copy(),
  varm = m2_50_PCA_test.varm.copy(),
  layers = m2_50_PCA_test.layers.copy(),
  raw = m2_50_PCA_test.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = m2_50_PCA_test.obsp.copy(),
  varp = m2_50_PCA_test.varp
  )
b_test.__dict__['_raw'].__dict__['_var'] = m2_50_PCA_test.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 2 data object
b_test.write(par_m2_test)

############################################################################################################################################
# Method 3:
############################################################################################################################################
## Copy the data object
m3_50_NMF_factors_train = train_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m3_50_NMF_factors_train)

## Normalize
m3_50_NMF_factors_train = normalize(m3_50_NMF_factors_train)

## Logarithmize
m3_50_NMF_factors_train = logarithmize(m3_50_NMF_factors_train)

## NMF model
model = NMF(n_components=n_rank, init='random', random_state=0, max_iter = 10000) 

## Fit the model
model.fit(m3_50_NMF_factors_train.X)

## W = transformed data matrix, V = original feature matrix
W = model.transform(m3_50_NMF_factors_train.X)
H = model.components_
W.shape
H.shape

H_df = pd.DataFrame(H, columns=m3_50_NMF_factors_train.var_names)
H_df

W_df = pd.DataFrame(W)
W_df

## Rename disease ontology
m3_50_NMF_factors_train.obs['disease__ontology_label'] = m3_50_NMF_factors_train.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## Check class
type(W), type(m3_50_NMF_factors_train.X)

## Create dataframe
new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])

## Create new anndata object for method 3
c_train = sc.AnnData(X = W.copy(),
  obs = m3_50_NMF_factors_train.obs.copy(),
  var = new_var,
  uns = m3_50_NMF_factors_train.uns.copy(),
  obsm = m3_50_NMF_factors_train.obsm.copy(),
  varm = m3_50_NMF_factors_train.varm.copy(),
  layers = m3_50_NMF_factors_train.layers.copy(),
  raw = m3_50_NMF_factors_train.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = m3_50_NMF_factors_train.obsp.copy(),
  varp = m3_50_NMF_factors_train.varp
  )
c_train.__dict__['_raw'].__dict__['_var'] = m3_50_NMF_factors_train.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 3 data object
c_train.write(par_m3_train)

## Check the reconstruction error -- it is okay ish
model.reconstruction_err_

top_genes_dict = {}
for i in range(H_df.shape[0]):  # Iterate through each factor
    top_genes = H_df.iloc[i].nlargest(H_df.shape[1])
    top_genes_dict[f'feature_{i}'] = top_genes
top_genes_df = pd.DataFrame(top_genes_dict)
top_genes_df.to_csv(par_top_genes_factor)


## Apply NMF model to test set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Copy the data object
m3_50_NMF_factors_test = test_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m3_50_NMF_factors_test)

## Normalize
m3_50_NMF_factors_test = normalize(m3_50_NMF_factors_test)

## Logarithmize
m3_50_NMF_factors_test = logarithmize(m3_50_NMF_factors_test)

## NMF model
W = model.transform(m3_50_NMF_factors_test.X)
H = model.components_
W.shape
H.shape

H_df = pd.DataFrame(H, columns=m3_50_NMF_factors_test.var_names)
H_df

W_df = pd.DataFrame(W)
W_df

## Rename disease ontology
m3_50_NMF_factors_test.obs['disease__ontology_label'] = m3_50_NMF_factors_test.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## Check class
type(W), type(m3_50_NMF_factors_test.X)

## Create dataframe
new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])

## Create new anndata object for method 3
c_test = sc.AnnData(X = W.copy(),
  obs = m3_50_NMF_factors_test.obs.copy(),
  var = new_var,
  uns = m3_50_NMF_factors_test.uns.copy(),
  obsm = m3_50_NMF_factors_test.obsm.copy(),
  varm = m3_50_NMF_factors_test.varm.copy(),
  layers = m3_50_NMF_factors_test.layers.copy(),
  raw = m3_50_NMF_factors_test.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = m3_50_NMF_factors_test.obsp.copy(),
  varp = m3_50_NMF_factors_test.varp
  )
c_test.__dict__['_raw'].__dict__['_var'] = m3_50_NMF_factors_test.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 3 data object
c_test.write(par_m3_test)


############################################################################################################################################
# Method 4: scETM 
############################################################################################################################################
## Train topic model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Copy the data object
m4_50_ETM_topics_train = train_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m4_50_ETM_topics_train)

## Set the necesarry column
m4_50_ETM_topics_train.obs["batch_indices"] = m4_50_ETM_topics_train.obs["orig.ident"]
m4_50_ETM_topics_train.obs["assigned_cluster"] = m4_50_ETM_topics_train.obs["Cell_Type"]

## scETM model 
model = scETM(m4_50_ETM_topics_train.n_vars, m4_50_ETM_topics_train.obs.batch_indices.nunique())

## train the model 
trainer = UnsupervisedTrainer(model, m4_50_ETM_topics_train, test_ratio=0.1)
trainer.train(n_epochs = 100, eval_every = 50, eval_kwargs = dict(cell_type_col = 'assigned_cluster'), save_model_ckpt = False) ## increase 12000

## Get cell embeddings
model.get_all_embeddings_and_nll(m4_50_ETM_topics_train)
delta, alpha, rho = map(pd.DataFrame, [m4_50_ETM_topics_train.obsm['delta'], m4_50_ETM_topics_train.uns['alpha'], m4_50_ETM_topics_train.varm['rho']])
delta.index = m4_50_ETM_topics_train.obs_names
rho.index = m4_50_ETM_topics_train.var_names
delta.shape, alpha.shape, rho.shape

## Get top 200 genes per topic
beta = rho @ alpha.T  # (gene, topic)
top_words = pd.DataFrame(m4_50_ETM_topics_train.var_names.values[np.argsort(beta.values, axis=0)[:-201:-1]])  # (n_top, topic)
top_words.to_csv(par_top_genes_topic)

## save AnnData object
d_train = sc.AnnData(X=m4_50_ETM_topics_train.X.copy(),
                  obs=m4_50_ETM_topics_train.obs.copy(),
                  var=m4_50_ETM_topics_train.var.copy(),
                  uns=m4_50_ETM_topics_train.uns.copy(),
                  obsm=m4_50_ETM_topics_train.obsm.copy(),
                  varm=m4_50_ETM_topics_train.varm.copy(),
                  layers=m4_50_ETM_topics_train.layers.copy(),
                  raw=m4_50_ETM_topics_train.raw.copy(),
                  dtype="float32",
                  shape=None,
                  obsp=m4_50_ETM_topics_train.obsp.copy(),
                  varp=m4_50_ETM_topics_train.varp
                  )
d_train.__dict__['_raw'].__dict__['_var'] = m4_50_ETM_topics_train.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 4 data object
d_train.write(par_m4_train)

## Apply topic model to test set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Copy the data object
m4_50_ETM_topics_test = test_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m4_50_ETM_topics_test)

## Set the necesarry column
m4_50_ETM_topics_test.obs["batch_indices"] = m4_50_ETM_topics_test.obs["orig.ident"]
m4_50_ETM_topics_test.obs["assigned_cluster"] = m4_50_ETM_topics_test.obs["Cell_Type"]

## Get cell embeddings
model.get_all_embeddings_and_nll(m4_50_ETM_topics_test)
delta, alpha, rho = map(pd.DataFrame, [m4_50_ETM_topics_test.obsm['delta'], m4_50_ETM_topics_test.uns['alpha'], m4_50_ETM_topics_test.varm['rho']])
delta.index = m4_50_ETM_topics_test.obs_names
rho.index = m4_50_ETM_topics_test.var_names
delta.shape, alpha.shape, rho.shape

## save AnnData object
d_test = sc.AnnData(X=m4_50_ETM_topics_test.X.copy(),
                  obs=m4_50_ETM_topics_test.obs.copy(),
                  var=m4_50_ETM_topics_test.var.copy(),
                  uns=m4_50_ETM_topics_test.uns.copy(),
                  obsm=m4_50_ETM_topics_test.obsm.copy(),
                  varm=m4_50_ETM_topics_test.varm.copy(),
                  layers=m4_50_ETM_topics_test.layers.copy(),
                  raw=m4_50_ETM_topics_test.raw.copy(),
                  dtype="float32",
                  shape=None,
                  obsp=m4_50_ETM_topics_test.obsp.copy(),
                  varp=m4_50_ETM_topics_test.varp
                  )
d_test.__dict__['_raw'].__dict__['_var'] = m4_50_ETM_topics_test.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})

## Save method 4 data object
d_test.write(par_m4_test)
