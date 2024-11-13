## This script was used to split the Kamath object. We use astrocytes as an example. 

salloc -A def-sfarhan --time=0-5 -c 1 --mem=300g

module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  


###################################
# Adjustable parameters
###################################
par_ann_data = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Kamath.h5ad"

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
# Cell Specific parameters -- astro
###################################
par_keep_cell_type = "astro"
par_test_anndata = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/test_kamath_astro.h5ad'
par_train_anndata = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/train_kamath_astro.h5ad'

###################################
# Load anndata
###################################
adata = sc.read_h5ad(par_ann_data)
adata.obs
adata = adata[adata.obs['Cell_Type'] == par_keep_cell_type]
adata.obs['Cell_Type']
print("genes: ", adata.var_names) 
print("cells: ", adata.obs_names) 
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

###################################
# create new_annData object without LBD cell samples
###################################
adata.obs = adata.obs.reset_index() 
adata.obs

adata.obs['Disease_Status']

adata.obs
adata.obs.shape

set(adata.obs['Disease_Status'])

## split data set into test and train
train, test = train_test_split(adata, test_size=0.2, random_state=5)

train.obs.shape
test.obs.shape

## save test objct ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_anndata = sc.AnnData(X = test.raw.X.toarray(),
  obs = test.obs.copy(),
  var = test.var.copy(),
  uns = test.uns.copy(),
  obsm = test.obsm.copy(),
  varm = test.varm.copy(),
  layers = test.layers.copy(),
  raw = test.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = test.obsp.copy(),
  varp = test.varp
  )

np.max(test_anndata.X)

## Necessary for the object to work.
test_anndata.__dict__['_raw'].__dict__['_var'] = test.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

## save the new annData object
test_anndata.write(par_test_anndata)

## save train object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
train_anndata = sc.AnnData(X = train.raw.X.toarray(),
  obs = train.obs.copy(),
  var = train.var.copy(),
  uns = train.uns.copy(),
  obsm = train.obsm.copy(),
  varm = train.varm.copy(),
  layers = train.layers.copy(),
  raw = train.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = train.obsp.copy(),
  varp = train.varp
  )

np.max(train_anndata.X)

## Necessary for the object to work.
train_anndata.__dict__['_raw'].__dict__['_var'] = train.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

## save the new annData object
train_anndata.write(par_train_anndata)

