## This script was used to identify HVGs in the Wang object. We use astrocytes as an example. 

salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

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
# Adjustable parameters cell type specific -- astro
######################################################################
par_test_anndata = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/test_wang_astro.h5ad'
par_train_anndata = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/train_wang_astro.h5ad'

par_m1_test = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/HVG_test_wang_astro.h5ad'
par_m1_train = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/HVG_train_wang_astro.h5ad'


###################################
# Load anndata
###################################
train_anndata = sc.read_h5ad(par_train_anndata)
test_anndata =  sc.read_h5ad(par_test_anndata)

##################################################################################################################################################################################################################
# Method 1: top HVGs in > 3 cells (raw)
##################################################################################################################################################################################################################
## identify highly variable genes on train dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
m1_5000_genes_train.obs['Disease_Status'] = m1_5000_genes_train.obs['Disease_Status'].replace({'PD': 1, 'control': 0})

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
m1_5000_genes_test.obs['Disease_Status'] = m1_5000_genes_test.obs['Disease_Status'].replace({'PD': 1, 'control': 0})

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

