############################################################
############################################################
############################################################
############################################################ NOTES
# This script can be used to develop NN classifiers to distinguish individual cells according to treatment or disease status 
# and apply LIME to the models to reveal the most important features for the classification decision.

# Specifically, this script will:
#1) split the AnnData object into 80% for training and 20% for testing;
#2) identify highly variable genes for input into the NN models;
#3) train the NN models via 5 fold cross-validation;
#4) apply the models to the test set;
#5) decode the models using LIME.

# This script is compatible with AnnData objects that contain two treatment/disease status categories.
# If your AnnData object contains > two treatment/disease status categories, please subset your AnnData object prior to the analysis.

# This script has been tested in Python v3.8.10

############################################################
############################################################
############################################################
############################################################ Analytical parameters
## Before proceeding with the analysis, be sure to define the following variables:

# Path to the working directory. The outputs from this script will be placed there.
par_out_dir = "/home/fiorini9/scratch/ML_github_test/"

# Path to AnnData Scanpy object.
par_ann_data = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5ad" 

# Meta data column name that contains the cell type annotations.
par_cell_type_column = "Cell_Type"

# Cell type to be retained for this analysis.
par_keep_cell_type = "astro"

# Meta data column name that contains the treatment/disease status information.
par_treatment_disease_column = "Disease_Status"

# Control group label described in the treatment/disease status column. 
par_control = "ctrl"

# Case group label described in the treatment/disease status column. 
par_case = "ALS"

# Meta data column name that contains the IDs of the individual subjects.
par_subject_ID = "Sample_ID"

# Activation function for NN.
par_NN_activation = "relu"

# Solver for NN.
par_NN_solver = "adam"

# Maximum interation for NN.
par_NN_max_iteration = 500


############################################################
############################################################
############################################################
############################################################ Load required libraries
# Standard library imports
import sys
import random
import time
from math import sqrt
from pprint import pprint

# Data manipulation and analysis
import numpy as np
import pandas as pd
import copy
import re

# ML libraries
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, cross_val_score, KFold, GridSearchCV, RandomizedSearchCV
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report, f1_score, 
                             balanced_accuracy_score, make_scorer, mean_squared_error, r2_score, 
                             ConfusionMatrixDisplay)
import sklearn.metrics as metrics
from sklearn.neural_network import MLPClassifier, MLPRegressor

# Other packages
import anndata as ad
import lime
import lime.lime_tabular
import importlib


############################################################
############################################################
############################################################
############################################################ Test-train split

###################################
# Preliminary parameters
###################################
## Set paths to save test/train objects
par_test_anndata = f"{par_out_dir}test_{par_keep_cell_type}_anndata.h5ad"
par_train_anndata = f"{par_out_dir}train_{par_keep_cell_type}_anndata.h5ad"

###################################
# Prepare input AnnData
###################################
## Load
adata = sc.read_h5ad(par_ann_data)
adata = adata[adata.obs[par_cell_type_column] == par_keep_cell_type]

## Check if the cell type subset worked 
unique_values = adata.obs[par_cell_type_column].unique()

## Check if all unique values are equal to the predefined value
if np.all(unique_values == par_keep_cell_type):
    print("Cell type subset completed successfully.")
else:
    raise ValueError("Cell type subset failed.")

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

###################################
# Create and save train/test AnnData objects
###################################
## reset index
adata.obs = adata.obs.reset_index() 

## split data set into test and train
train, test = train_test_split(adata, test_size=0.2, random_state=5)

## Prepare test object
test_anndata = sc.AnnData(X = test.X,
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

test_anndata.__dict__['_raw'].__dict__['_var'] = test.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

## Save test object
test_anndata.write(par_test_anndata)

## Prepare train object
train_anndata = sc.AnnData(X = train.X,
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

train_anndata.__dict__['_raw'].__dict__['_var'] = train.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

## Save the train object
train_anndata.write(par_train_anndata)


############################################################
############################################################
############################################################
############################################################ HVG feature selection

###################################
# Define helper functions
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

###################################
# Preliminary parameters
###################################
## Previously saved test/train anndata objects
par_test_anndata = f"{par_out_dir}test_{par_keep_cell_type}_anndata.h5ad"
par_train_anndata = f"{par_out_dir}train_{par_keep_cell_type}_anndata.h5ad"

## Set paths to save test/train objects reduced to HVGs
par_m1_test = f"{par_out_dir}HVG_test_{par_keep_cell_type}_anndata.h5ad"
par_m1_train = f"{par_out_dir}HVG_train_{par_keep_cell_type}_anndata.h5ad"

## Set paths to save test/train objects metadata for downstream processing
par_test_meta = f"{par_out_dir}test_metadata_{par_keep_cell_type}.csv"
par_train_meta = f"{par_out_dir}train_metadata_{par_keep_cell_type}.csv"

###################################
# Load test/train AnnData objects
###################################
## Load previously save AnnData Objects
train_anndata = sc.read_h5ad(par_train_anndata)
test_anndata =  sc.read_h5ad(par_test_anndata)

###################################
# Identify highly variable genes on train dataset
###################################
## Copy the trains object
m1_5000_genes_train = train_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m1_5000_genes_train)

## Normalize data
m1_5000_genes_train = normalize(m1_5000_genes_train)

## Logarithmize data
m1_5000_genes_train = logarithmize(m1_5000_genes_train)

## Identify highly variable genes
m1_5000_genes_train = HVGs(m1_5000_genes_train)

## Getting HVGs and store in dataframe
print('size before HVGs: ' + str(m1_5000_genes_train.X.shape) + str(m1_5000_genes_train.obs.shape))
m1_5000_genes_train = select_HVGs(m1_5000_genes_train)
print('size after HVGs: ' + str(m1_5000_genes_train.X.shape) + str(m1_5000_genes_train.obs.shape))

## Rename disease ontology
m1_5000_genes_train.obs[par_treatment_disease_column] = m1_5000_genes_train.obs[par_treatment_disease_column].replace({par_case: 1, par_control: 0})

## Prepare train object
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

## Save train object
a_train.write(par_m1_train)

## Print train set metdata to be used for LIME processing downstream      
a_train.obs.to_csv(par_train_meta)


###################################
# pull out the same highly variable genes on test dataset
###################################
## copy the data object
m1_5000_genes_test = test_anndata.copy()

## Annotate the MT genes
annotate_mt_genes(m1_5000_genes_test)

## Normalize data
m1_5000_genes_test = normalize(m1_5000_genes_test)

## Logarithmize data
m1_5000_genes_test = logarithmize(m1_5000_genes_test)

## only keep HGV gnes from train datasets
m1_5000_genes_test = m1_5000_genes_test[:,m1_5000_genes_train.var_names] ## select same HVGs from train

## Rename disease ontology
m1_5000_genes_test.obs[par_treatment_disease_column] = m1_5000_genes_test.obs[par_treatment_disease_column].replace({par_case: 1, par_control: 0}) ## adjust disease status meta data column (i.e Disease_Status) as well as the actual values (i.e PD or Control)

## Prepare test object
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

## Save test object
a_test.write(par_m1_test)

## Print test set metdata to be used for LIME processing downstream      
a_test.obs.to_csv(par_test_meta)


############################################################
############################################################
############################################################
############################################################ Train NN models

###################################
# Preliminary parameters
###################################
## Previously saved test/train anndata objects
par_m1_test = f"{par_out_dir}HVG_test_{par_keep_cell_type}_anndata.h5ad"
par_m1_train = f"{par_out_dir}HVG_train_{par_keep_cell_type}_anndata.h5ad"

## Path to dataframe describing the performance metrics on the test set
par_performance = f"{par_out_dir}test_performance_{par_keep_cell_type}.csv"
par_CM = f"{par_out_dir}test_confusion_matrix_{par_keep_cell_type}.csv"

## Path to dataframe describing the cell indices of the test set - important for downstream LIME computations
NN_par_cell_index_fold1 = f"{par_out_dir}cell_index_fold1_{par_keep_cell_type}.csv"
NN_par_cell_index_fold2 = f"{par_out_dir}cell_index_fold2_{par_keep_cell_type}.csv"
NN_par_cell_index_fold3 = f"{par_out_dir}cell_index_fold3_{par_keep_cell_type}.csv"
NN_par_cell_index_fold4 = f"{par_out_dir}cell_index_fold4_{par_keep_cell_type}.csv"
NN_par_cell_index_fold5 = f"{par_out_dir}cell_index_fold5_{par_keep_cell_type}.csv"

## Path to dataframe containing LIME outputs for each of the 5-fold cross validation models on test set
NN_par_LIME_out_m1_fold1 = f"{par_out_dir}LIME_fold1_{par_keep_cell_type}.csv"
NN_par_LIME_out_m1_fold2 = f"{par_out_dir}LIME_fold2_{par_keep_cell_type}.csv"
NN_par_LIME_out_m1_fold3 = f"{par_out_dir}LIME_fold3_{par_keep_cell_type}.csv"
NN_par_LIME_out_m1_fold4 = f"{par_out_dir}LIME_fold4_{par_keep_cell_type}.csv"
NN_par_LIME_out_m1_fold5 = f"{par_out_dir}LIME_fold5_{par_keep_cell_type}.csv"

###################################
# preprocessing adata objects
###################################
## Load train data object
m1 = sc.read_h5ad(par_m1_train)

###################################
# 5-fold cross validation 
###################################
## Define n splits
kf = KFold(n_splits=5)

## Set variables
m1_X = m1.X
m1_y = m1.obs[par_treatment_disease_column]
kf.get_n_splits(m1_X)

## Fold 1
xx = list(kf.split(m1_X))
xx0 = list(xx[0])

train_index_fold1 = xx0[0]
test_index_fold1 = xx0[1]

m1_X_train_fold1, m1_X_test_fold1, m1_y_train_fold1, m1_y_test_fold1 = m1_X[train_index_fold1], m1_X[test_index_fold1], m1_y[train_index_fold1], m1_y[test_index_fold1]

## Fold 2
xx = list(kf.split(m1_X))
xx0 = list(xx[1])

train_index_fold2 = xx0[0]
test_index_fold2 = xx0[1]

m1_X_train_fold2, m1_X_test_fold2, m1_y_train_fold2, m1_y_test_fold2 = m1_X[train_index_fold2], m1_X[test_index_fold2], m1_y[train_index_fold2], m1_y[test_index_fold2]

## Fold 3
xx = list(kf.split(m1_X))
xx0 = list(xx[2])

train_index_fold3 = xx0[0]
test_index_fold3 = xx0[1]

m1_X_train_fold3, m1_X_test_fold3, m1_y_train_fold3, m1_y_test_fold3 = m1_X[train_index_fold3], m1_X[test_index_fold3], m1_y[train_index_fold3], m1_y[test_index_fold3]

## Fold 4
xx = list(kf.split(m1_X))
xx0 = list(xx[3])

train_index_fold4 = xx0[0]
test_index_fold4 = xx0[1]

m1_X_train_fold4, m1_X_test_fold4, m1_y_train_fold4, m1_y_test_fold4 = m1_X[train_index_fold4], m1_X[test_index_fold4], m1_y[train_index_fold4], m1_y[test_index_fold4]

## Fold 5
xx = list(kf.split(m1_X))
xx0 = list(xx[4])

train_index_fold5 = xx0[0]
test_index_fold5 = xx0[1]

m1_X_train_fold5, m1_X_test_fold5, m1_y_train_fold5, m1_y_test_fold5 = m1_X[train_index_fold5], m1_X[test_index_fold5], m1_y[train_index_fold5], m1_y[test_index_fold5]


###################################
# Define helper functions
###################################
def fit(X_train, y_train):
    model.fit(X_train, y_train)
    return model

def validation(model, X_val, y_val):
    prediction_val = model.predict(X_val)
    accuracy_val = metrics.accuracy_score(y_val, prediction_val)
    print("Validation accuracy = ", metrics.accuracy_score(y_val, prediction_val))
    return prediction_val, accuracy_val

def predict(model, X_test, y_test):
    prediction_test = model.predict(X_test)
    accuracy_test = metrics.accuracy_score(y_test, prediction_test)
    print("Test accuracy = ", metrics.accuracy_score(y_test, prediction_test))
    return metrics.accuracy_score(y_test, prediction_test), prediction_test

###################################
# Produce NN model and test on validation -- 5-fold cross validation
###################################
## Fold 1
MLP_m1_fold1 = MLPClassifier( activation=par_NN_activation, solver=par_NN_solver, max_iter=par_NN_max_iteration)
MLP_m1_fold1.fit(m1_X_train_fold1, m1_y_train_fold1)
MLP_m1_classifier_fold1 = MLP_m1_fold1.fit(m1_X_train_fold1, m1_y_train_fold1)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold1, m1_X_test_fold1, m1_y_test_fold1) 
## Fold 2
MLP_m1_fold2 = MLPClassifier( activation=par_NN_activation, solver=par_NN_solver, max_iter=par_NN_max_iteration)
MLP_m1_fold2.fit(m1_X_train_fold2, m1_y_train_fold2)
MLP_m1_classifier_fold2 = MLP_m1_fold2.fit(m1_X_train_fold2, m1_y_train_fold2)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold2, m1_X_test_fold2, m1_y_test_fold2) 
## Fold 3
MLP_m1_fold3 = MLPClassifier( activation=par_NN_activation, solver=par_NN_solver, max_iter=par_NN_max_iteration)
MLP_m1_fold3.fit(m1_X_train_fold3, m1_y_train_fold3)
MLP_m1_classifier_fold3 = MLP_m1_fold3.fit(m1_X_train_fold3, m1_y_train_fold3)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold3, m1_X_test_fold3, m1_y_test_fold3) 
## Fold 4
MLP_m1_fold4 = MLPClassifier( activation=par_NN_activation, solver=par_NN_solver, max_iter=par_NN_max_iteration)
MLP_m1_fold4.fit(m1_X_train_fold4, m1_y_train_fold4)
MLP_m1_classifier_fold4 = MLP_m1_fold4.fit(m1_X_train_fold4, m1_y_train_fold4)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold4, m1_X_test_fold4, m1_y_test_fold4) 
## Fold 5
MLP_m1_fold5 = MLPClassifier( activation=par_NN_activation, solver=par_NN_solver, max_iter=par_NN_max_iteration)
MLP_m1_fold5.fit(m1_X_train_fold5, m1_y_train_fold5)
MLP_m1_classifier_fold5 = MLP_m1_fold5.fit(m1_X_train_fold5, m1_y_train_fold5)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold5, m1_X_test_fold5, m1_y_test_fold5) 

############################################################
############################################################
############################################################
############################################################ Model evaluation on test set

###################################
# load test objects
###################################
## Load test data object
m1_test = sc.read_h5ad(par_m1_test)

## Define the variables
X_test_m1 = m1_test.X
y_test_m1 = m1_test.obs[par_treatment_disease_column]

## Retrieve number of HVGs
num_genes = X_test_m1.shape[1]


###################################
# Apply model to test dataset -- Fold 1
###################################
## Prediction
prediction_test_fold1 = MLP_m1_fold1.predict(X_test_m1)
predictions = MLP_m1_classifier_fold1.predict_proba(X_test_m1)

## Accuracy
accuracy_test_m1_fold1 = metrics.accuracy_score(y_test_m1, prediction_test_fold1)

## Balanced accuracy
BA_test_m1_fold1 = balanced_accuracy_score(y_test_m1, prediction_test_fold1)

## F1
F1_test_m1_fold1 = f1_score(y_test_m1, prediction_test_fold1, pos_label=1)

## Confusuon matrix 
cm_array = confusion_matrix(y_test_m1, prediction_test_fold1)
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': [par_case, par_control,par_case, par_control],
'predict_label': [par_case, par_control, par_control, par_case ],
'method': ['fold1','fold1','fold1','fold1']
}
confusion_matrix_m1_fold1 = pd.DataFrame(data) 

###################################
# Apply model to test dataset -- Fold 2
###################################
## Prediction
prediction_test_fold2 = MLP_m1_fold2.predict(X_test_m1)
predictions = MLP_m1_classifier_fold2.predict_proba(X_test_m1)

## Accuracy
accuracy_test_m1_fold2 = metrics.accuracy_score(y_test_m1, prediction_test_fold2)

## Balanced accuracy
BA_test_m1_fold2 = balanced_accuracy_score(y_test_m1, prediction_test_fold2)

## F1
F1_test_m1_fold2 = f1_score(y_test_m1, prediction_test_fold2, pos_label=1)

## Confusuon matrix 
cm_array = confusion_matrix(y_test_m1, prediction_test_fold2)
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': [par_case, par_control,par_case, par_control],
'predict_label': [par_case, par_control, par_control, par_case ],
'method': ['fold2','fold2','fold2','fold2']
}
confusion_matrix_m1_fold2 = pd.DataFrame(data) 


###################################
# Apply model to test dataset -- Fold 3
###################################
## Prediction
prediction_test_fold3 = MLP_m1_fold3.predict(X_test_m1)
predictions = MLP_m1_classifier_fold3.predict_proba(X_test_m1)

## Accuracy
accuracy_test_m1_fold3 = metrics.accuracy_score(y_test_m1, prediction_test_fold3)

## Balanced accuracy
BA_test_m1_fold3 = balanced_accuracy_score(y_test_m1, prediction_test_fold3)

## F1
F1_test_m1_fold3 = f1_score(y_test_m1, prediction_test_fold3, pos_label=1)

## Confusuon matrix 
cm_array = confusion_matrix(y_test_m1, prediction_test_fold3)
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': [par_case, par_control,par_case, par_control],
'predict_label': [par_case, par_control, par_control, par_case ],
'method': ['fold3','fold3','fold3','fold3']
}
confusion_matrix_m1_fold3 = pd.DataFrame(data) 

###################################
# Apply model to test dataset -- Fold 4
###################################
## Prediction
prediction_test_fold4 = MLP_m1_fold4.predict(X_test_m1)
predictions = MLP_m1_classifier_fold4.predict_proba(X_test_m1)

## Accuracy
accuracy_test_m1_fold4 = metrics.accuracy_score(y_test_m1, prediction_test_fold4)

## Balanced accuracy
BA_test_m1_fold4 = balanced_accuracy_score(y_test_m1, prediction_test_fold4)

## F1
F1_test_m1_fold4 = f1_score(y_test_m1, prediction_test_fold4, pos_label=1)

## Confusuon matrix 
cm_array = confusion_matrix(y_test_m1, prediction_test_fold4)
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': [par_case, par_control,par_case, par_control],
'predict_label': [par_case, par_control, par_control, par_case ],
'method': ['fold4','fold4','fold4','fold4']
}
confusion_matrix_m1_fold4 = pd.DataFrame(data) 

###################################
# Apply model to test dataset -- Fold 5
###################################
## Prediction
prediction_test_fold5 = MLP_m1_fold5.predict(X_test_m1)
predictions = MLP_m1_classifier_fold5.predict_proba(X_test_m1)

## Accuracy
accuracy_test_m1_fold5 = metrics.accuracy_score(y_test_m1, prediction_test_fold5)

## Balanced accuracy
BA_test_m1_fold5 = balanced_accuracy_score(y_test_m1, prediction_test_fold5)

## F1
F1_test_m1_fold5 = f1_score(y_test_m1, prediction_test_fold5, pos_label=1)

## Confusuon matrix 
cm_array = confusion_matrix(y_test_m1, prediction_test_fold5)
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': [par_case, par_control,par_case, par_control],
'predict_label': [par_case, par_control, par_control, par_case ],
'method': ['fold4','fold4','fold4','fold4']
}
confusion_matrix_m1_fold5 = pd.DataFrame(data) 


####################################
## Test set performance metrics and confusion matrix
####################################
## Assemble performace metrics dataframe
performance = {'model': ["fold1", "fold2","fold3","fold4","fold5",],
'accuracy': [accuracy_test_m1_fold1, accuracy_test_m1_fold2, accuracy_test_m1_fold3, accuracy_test_m1_fold4, accuracy_test_m1_fold5],
'balanced_accuracy': [BA_test_m1_fold1, BA_test_m1_fold2, BA_test_m1_fold3, BA_test_m1_fold4, BA_test_m1_fold5],
'F1': [F1_test_m1_fold1, F1_test_m1_fold2, F1_test_m1_fold3, F1_test_m1_fold4, F1_test_m1_fold5]
}

## Print performace metrics dataframe
performance_df = pd.DataFrame(performance) 
performance_df.to_csv(par_performance)

## Assemble confusion matrix dataframe
all_confusion_matrices = pd.concat([confusion_matrix_m1_fold1, 
                                    confusion_matrix_m1_fold2, 
                                    confusion_matrix_m1_fold3, 
                                    confusion_matrix_m1_fold4, 
                                    confusion_matrix_m1_fold5], 
                                   ignore_index=True)

## Print confusion matrix dataframe
all_confusion_matrices.to_csv(par_CM)


####################################
## Print the cell indices of the test set 
####################################
## Fold 1
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test_fold1
df_index['cell_index'] = df_index.index
df_index.to_csv(NN_par_cell_index_fold1)

## Fold 2
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test_fold2
df_index['cell_index'] = df_index.index
df_index.to_csv(NN_par_cell_index_fold2)

## Fold 3
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test_fold3
df_index['cell_index'] = df_index.index
df_index.to_csv(NN_par_cell_index_fold3)

## Fold 4
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test_fold4
df_index['cell_index'] = df_index.index
df_index.to_csv(NN_par_cell_index_fold4)

## Fold 5
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test_fold5
df_index['cell_index'] = df_index.index
df_index.to_csv(NN_par_cell_index_fold5)


####################################
## Prepare data for LIME
####################################
## Convert to dense array
X_test_dense = X_test_m1.toarray() if hasattr(X_test_m1, 'toarray') else X_test_m1

####################################
## Apply LIME to NN - fold 1
####################################
## Convert to dense array
m1_X_train_fold1_dense = m1_X_train_fold1.toarray() if hasattr(m1_X_train_fold1, 'toarray') else m1_X_train_fold1

## Convert array to dataframe
df=pd.DataFrame(m1_X_train_fold1_dense,) 
df_test=pd.DataFrame(X_test_dense,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## Initiate Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=[par_control, par_case], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## LIME loop
df = pd.DataFrame()
df[0] = []
df[1] = []
df['cell_index'] = []

for i in range(0, num_cells):
    exp = explainer.explain_instance(df_test.values[i], MLP_m1_classifier_fold1.predict_proba, num_features=num_genes)
    dd = exp.as_list()
    df2=pd.DataFrame(dd,) 
    df2['cell_index'] = df_test.iloc[[i]].index[0]
    df=pd.concat([df, df2])

## Save LIME outputs
df.to_csv(NN_par_LIME_out_m1_fold1)

####################################
## Apply LIME to NN - fold 2
####################################
## Convert to dense array
m1_X_train_fold2_dense = m1_X_train_fold2.toarray() if hasattr(m1_X_train_fold2, 'toarray') else m1_X_train_fold2

## Convert array to dataframe
df=pd.DataFrame(m1_X_train_fold2_dense,) 
df_test=pd.DataFrame(X_test_dense,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## Initiate Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=[par_control, par_case], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## LIME loop
df = pd.DataFrame()
df[0] = []
df[1] = []
df['cell_index'] = []

for i in range(0, num_cells):
    exp = explainer.explain_instance(df_test.values[i], MLP_m1_classifier_fold2.predict_proba, num_features=num_genes)
    dd = exp.as_list()
    df2=pd.DataFrame(dd,) 
    df2['cell_index'] = df_test.iloc[[i]].index[0]
    df=pd.concat([df, df2])

## Save LIME outputs
df.to_csv(NN_par_LIME_out_m1_fold2)

####################################
## Apply LIME to NN - fold 3
####################################
## Convert to dense array
m1_X_train_fold3_dense = m1_X_train_fold3.toarray() if hasattr(m1_X_train_fold3, 'toarray') else m1_X_train_fold3

## Convert array to dataframe
df=pd.DataFrame(m1_X_train_fold3_dense,) 
df_test=pd.DataFrame(X_test_dense,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## Initiate Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=[par_control, par_case], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## LIME loop
df = pd.DataFrame()
df[0] = []
df[1] = []
df['cell_index'] = []

for i in range(0, num_cells):
    exp = explainer.explain_instance(df_test.values[i], MLP_m1_classifier_fold3.predict_proba, num_features=num_genes)
    dd = exp.as_list()
    df2=pd.DataFrame(dd,) 
    df2['cell_index'] = df_test.iloc[[i]].index[0]
    df=pd.concat([df, df2])

## Save LIME outputs
df.to_csv(NN_par_LIME_out_m1_fold3)

####################################
## Apply LIME to NN - fold 4
####################################
## Convert to dense array
m1_X_train_fold4_dense = m1_X_train_fold4.toarray() if hasattr(m1_X_train_fold4, 'toarray') else m1_X_train_fold4

## Convert array to dataframe
df=pd.DataFrame(m1_X_train_fold4_dense,) 
df_test=pd.DataFrame(X_test_dense,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## Initiate Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['ctrl', 'ALS'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## LIME loop
df = pd.DataFrame()
df[0] = []
df[1] = []
df['cell_index'] = []

for i in range(0, num_cells):
    exp = explainer.explain_instance(df_test.values[i], MLP_m1_classifier_fold4.predict_proba, num_features=num_genes)
    dd = exp.as_list()
    df2=pd.DataFrame(dd,) 
    df2['cell_index'] = df_test.iloc[[i]].index[0]
    df=pd.concat([df, df2])

## Save LIME outputs
df.to_csv(NN_par_LIME_out_m1_fold4)

####################################
## Apply LIME to NN - fold 5
####################################
## Convert to dense array
m1_X_train_fold5_dense = m1_X_train_fold5.toarray() if hasattr(m1_X_train_fold5, 'toarray') else m1_X_train_fold5

## Convert array to dataframe
df=pd.DataFrame(m1_X_train_fold5_dense,) 
df_test=pd.DataFrame(X_test_dense,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## Initiate Explainer
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['ctrl', 'ALS'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## LIME loop
df = pd.DataFrame()
df[0] = []
df[1] = []
df['cell_index'] = []

for i in range(0, num_cells):
    exp = explainer.explain_instance(df_test.values[i], MLP_m1_classifier_fold5.predict_proba, num_features=num_genes)
    dd = exp.as_list()
    df2=pd.DataFrame(dd,) 
    df2['cell_index'] = df_test.iloc[[i]].index[0]
    df=pd.concat([df, df2])

## Save LIME outputs
df.to_csv(NN_par_LIME_out_m1_fold5)



############################################################
############################################################
############################################################
############################################################ LIME outputs processing
###################################
## Fold evaluation and selection
###################################
## Read in performance
performance_df = pd.read_csv(par_performance)

## Find the top performing model
top_fold_row = performance_df.loc[performance_df['balanced_accuracy'].idxmax()]
top_fold = top_fold_row['model']


###################################
## Load AnnData and select cell type
###################################
## Load
adata = sc.read_h5ad(par_ann_data)
adata = adata[adata.obs[par_cell_type_column] == par_keep_cell_type]

## Check if the cell type subset worked 
unique_values = adata.obs[par_cell_type_column].unique()

## Check if all unique values are equal to the predefined value
if np.all(unique_values == par_keep_cell_type):
    print("Cell type subset completed successfully.")
else:
    raise ValueError("Cell type subset failed.")

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

## reset index
adata.obs = adata.obs.reset_index() 


###################################
## Load metadata to perform sample-based processing
###################################
## Read in test object metadata
meta_data = pd.read_csv(par_test_meta)

## rename metadata
test = meta_data


####################################
## LIME processing
###################################
## Construct paths to relevant files for top fold
NN_LIME_fold = f"{par_out_dir}LIME_{top_fold}_{par_keep_cell_type}.csv"
NN_cell_index_fold = f"{par_out_dir}cell_index_{top_fold}_{par_keep_cell_type}.csv"

## Read in LIME output
LIME = pd.read_csv(NN_LIME_fold)

## Convert column names
LIME.columns = ["X", "X0", "importance", "cell_index"]

## Process feature importance values
LIME['test'] = LIME['X0'].str.count(' ')
LIME.loc[LIME['test'] == 4, 'feature'] = LIME.loc[LIME['test'] == 4, 'X0'].str.extract(r'.* .* (.*) .* .*', expand=False)
LIME.loc[LIME['test'] == 2, 'feature'] = LIME.loc[LIME['test'] == 2, 'X0'].str.extract(r'(.*) .* .*', expand=False)

## Merge with cell index and only retain correctly classified cells
index = pd.read_csv(NN_cell_index_fold)
index = index[[par_treatment_disease_column, 'predicted_label', 'cell_index']]
LIME = LIME[['X0', 'importance', 'cell_index', 'feature']]
old_length = LIME.shape[0]
LIME = pd.merge(LIME, index, on='cell_index', how='outer')
new_length = LIME.shape[0]

## Check to ensure merge was successful
if old_length != new_length:
    raise ValueError("The merge was unsuccessful.")
else:
    print("The merge was successful.")

## Only keep correctly classified instances
LIME = LIME[LIME[par_treatment_disease_column] == LIME['predicted_label']]

## Subset the metadata object -- remove incorrectly classified cells from test set
test = test.rename(columns={'Unnamed: 0': 'X'})
test_2 = test[test['X'].isin(LIME['cell_index'])]
test_2 = test_2[['X', par_subject_ID]]

if len(test_2) == len(LIME['cell_index'].unique()):
    print("Row count matches unique cell_index count before duplicates removed.")
else:
    print("Row count does not match unique cell_index count before duplicates removed.")

## Merge LIME dataframe
LIME_merge = pd.merge(LIME, test_2, left_on='cell_index', right_on='X', how='left')

if len(LIME_merge) == len(LIME):
    print("The number of rows in LIME_merge is the same as LIME.")
else:
    print("The number of rows in LIME_merge is not the same as LIME.")

## Sanity check to ensure Subject IDs were not disrupted during merge process
cross_tab = pd.crosstab(LIME_merge[par_treatment_disease_column], LIME_merge[par_subject_ID])
print(cross_tab)

unique_sample_ids = len(LIME_merge[par_subject_ID].unique())
print(f"Number of unique subjects: {unique_sample_ids}")

####################################
## Compute LIME feature importance Z-score
###################################
## Only retain treated/case cells
LIME_merge_PD = LIME_merge[LIME_merge[par_treatment_disease_column] == 1]

## Loop through each gene to compute importance for each individual subject
filler = []

for i in LIME_merge_PD[par_subject_ID].unique():
    temper = LIME_merge_PD[LIME_merge_PD[par_subject_ID] == i]
    LIME_avg = temper.groupby('feature', as_index=False).agg(
        Mean_feature_importance_donor=('importance', 'mean')
    )
    LIME_avg[par_subject_ID] = i
    data = np.abs(LIME_avg['Mean_feature_importance_donor'])
    z_scores_donor = (data - data.mean()) / data.std()
    LIME_avg['z_scores_donor'] = z_scores_donor
    filler.append(LIME_avg)

filler = pd.concat(filler, ignore_index=True)

## Sanity check to ensure that gene Z-score for individual subjects completed successfully
unique_feature_count = filler['feature'].nunique()
unique_Subject_ID = filler['Sample_ID'].nunique()

if unique_feature_count*unique_Subject_ID == len(filler):
    print("Subject-specific Z-score computed sucessfully.")
else:
    print("Subject-specific Z-score not computed sucessfully.")

## Compute average Z-score per gene across all subjects
filler_1 = filler.groupby('feature', as_index=False).agg(
    Mean_z_score_across_donors=('z_scores_donor', 'mean')
)
filler_1 = filler_1.sort_values(by='Mean_z_score_across_donors', ascending=False)
Complete_df = pd.merge(filler, filler_1, on='feature', how='left')
Complete_df = Complete_df.sort_values(by='Mean_z_score_across_donors', ascending=False)

## Sanity check to ensure that gene Z-score across all subjects completed successfully
if len(Complete_df) == len(filler):
    print("Gene Z-score across subjects computed sucessfully.")
else:
    print("Gene Z-score across subjects not computed sucessfully.")

## Save file
Complete_df.to_csv(f"{par_out_dir}LIME_unweighted_{par_keep_cell_type}.csv", index=False)

####################################
## Compute LIME feature importance Z-score weighted by percent expression
###################################
## Only retain treated/case cells
LIME_merge_PD = LIME_merge[LIME_merge[par_treatment_disease_column] == 1]

## Compute percent expression
express_fill_list = []
unique_features = LIME_merge_PD['feature'].unique()

for feature in unique_features:
    expression = adata.X[:, adata.var_names == feature].toarray().flatten()
    cells_expressing = expression[expression > 0]
    percent_expressed = len(cells_expressing) / adata.shape[0]
    express_fill_temp = pd.DataFrame({'percent_expressed': [percent_expressed], 'feature': [feature]})
    express_fill_list.append(express_fill_temp)

express_fill = pd.concat(express_fill_list, ignore_index=True)
express_fill['normalized_percent_expressed'] = express_fill['percent_expressed'] / express_fill['percent_expressed'].sum()

## Loop through each gene to compute importance for each individual subject
results_list = []

sample_id = '191112_ALS_112_snRNA-D1.RDS'
for sample_id in LIME_merge_PD[par_subject_ID].unique():
    temper = LIME_merge_PD[LIME_merge_PD[par_subject_ID] == sample_id]
    LIME_avg = temper.groupby('feature').agg(
        Mean_feature_importance_donor=('importance', 'mean')
    ).reset_index()
    LIME_avg[par_subject_ID] = sample_id
    LIME_avg = LIME_avg.merge(express_fill[['feature', 'normalized_percent_expressed']], on='feature', how='left')
    LIME_avg['Mean_feature_importance_donor_weighted'] = LIME_avg['Mean_feature_importance_donor'] * LIME_avg['normalized_percent_expressed']
    weighted_data = np.abs(LIME_avg['Mean_feature_importance_donor_weighted'])
    z_scores = (weighted_data - weighted_data.mean()) / weighted_data.std()
    LIME_avg['z_scores_weighted_donor'] = z_scores
    results_list.append(LIME_avg)

filler = pd.concat(results_list, ignore_index=True)

## Sanity check to ensure that gene Z-score for individual subjects completed successfully
unique_feature_count = filler['feature'].nunique()
unique_Subject_ID = filler['Sample_ID'].nunique()

if unique_feature_count*unique_Subject_ID == len(filler):
    print("Subject-specific Z-score computed sucessfully.")
else:
    print("Subject-specific Z-score not computed sucessfully.")

## Compute average Z-score per gene across subjects
filler_1 = filler.groupby('feature', as_index=False).agg(
    Mean_z_score_weighted_across_donor=('z_scores_weighted_donor', 'mean')
)
filler_1 = filler_1.sort_values(by='Mean_z_score_weighted_across_donor', ascending=False)
Complete_df_weighted = pd.merge(filler, filler_1, on='feature', how='left')
Complete_df_weighted = Complete_df_weighted.sort_values(by='Mean_z_score_weighted_across_donor', ascending=False)

## Sanity check to ensure that gene Z-score across all subjects completed successfully
if len(Complete_df_weighted) == len(filler):
    print("Gene Z-score across subjects computed sucessfully.")
else:
    print("Gene Z-score across subjects not computed sucessfully.")

## Save file
Complete_df_weighted.to_csv(f"{par_out_dir}LIME_weighted_{par_keep_cell_type}.csv", index=False)


