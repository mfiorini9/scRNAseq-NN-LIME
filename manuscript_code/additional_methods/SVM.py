## This script was used for logistic regression. 

salloc -A def-sfarhan --time=0-5 -c 1 --mem=100g

module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

###################################
# Adjustable parameters
###################################
## all methods
par_m1_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m1_HVGS_train.h5ad'
par_m2_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m2_PCA_train.h5ad' 
par_m3_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m3_NMF_factors_train.h5ad'
par_m4_train = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m4_ETM_train.h5ad'

par_m1_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m1_HVGS_test.h5ad'
par_m2_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m2_PCA_test.h5ad' 
par_m3_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m3_NMF_factors_test.h5ad'
par_m4_test = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/m4_ETM_test.h5ad'

## Support vector machine ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par_metric_summary_frame_SVM = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/SVM_5_fold_summary_eval.csv'
par_confusion_matrix_combo_SVM = '/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astro/SVM_confusion_matrix_summary_eval.csv'


###################################
# load library
###################################
## SVM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

import sys
import random
from sklearn import preprocessing 
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import importlib
required_libraries = ['torch', 'PIL', 'matplotlib', 
                      'numpy', 'pandas']
for lib in required_libraries:
    if importlib.util.find_spec(lib) is None:
        print("%s unavailable" % lib)

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset
import seaborn as sns
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn import metrics
import sklearn.metrics as metrics
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn import svm
import copy

from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, f1_score, balanced_accuracy_score, make_scorer
from sklearn.model_selection import cross_val_score

###################################
# preprocessing adata objects
###################################
## Load train data object
m1 = sc.read_h5ad(par_m1_train)
m2 = sc.read_h5ad(par_m2_train)
m3 = sc.read_h5ad(par_m3_train)
m4 = sc.read_h5ad(par_m4_train)

## fix disease ontology label for m4
m4.obs['disease__ontology_label'] = m4.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

###################################
# 5-fold cross validation for each method
###################################
kf = KFold(n_splits=5)

############################## Method 1 ##############################
m1_X = m1.X
m1_y = m1.obs['disease__ontology_label']
kf.get_n_splits(m1_X)

## fold 1
xx = list(kf.split(m1_X))
xx0 = list(xx[0])

train_index_fold1 = xx0[0]
test_index_fold1 = xx0[1]

m1_X_train_fold1, m1_X_test_fold1, m1_y_train_fold1, m1_y_test_fold1 = m1_X[train_index_fold1], m1_X[test_index_fold1], m1_y[train_index_fold1], m1_y[test_index_fold1]

## fold 2
xx = list(kf.split(m1_X))
xx0 = list(xx[1])

train_index_fold2 = xx0[0]
test_index_fold2 = xx0[1]

m1_X_train_fold2, m1_X_test_fold2, m1_y_train_fold2, m1_y_test_fold2 = m1_X[train_index_fold2], m1_X[test_index_fold2], m1_y[train_index_fold2], m1_y[test_index_fold2]

## fold 3
xx = list(kf.split(m1_X))
xx0 = list(xx[2])

train_index_fold3 = xx0[0]
test_index_fold3 = xx0[1]

m1_X_train_fold3, m1_X_test_fold3, m1_y_train_fold3, m1_y_test_fold3 = m1_X[train_index_fold3], m1_X[test_index_fold3], m1_y[train_index_fold3], m1_y[test_index_fold3]

## fold 4
xx = list(kf.split(m1_X))
xx0 = list(xx[3])

train_index_fold4 = xx0[0]
test_index_fold4 = xx0[1]

m1_X_train_fold4, m1_X_test_fold4, m1_y_train_fold4, m1_y_test_fold4 = m1_X[train_index_fold4], m1_X[test_index_fold4], m1_y[train_index_fold4], m1_y[test_index_fold4]

## fold 5
xx = list(kf.split(m1_X))
xx0 = list(xx[4])

train_index_fold5 = xx0[0]
test_index_fold5 = xx0[1]

m1_X_train_fold5, m1_X_test_fold5, m1_y_train_fold5, m1_y_test_fold5 = m1_X[train_index_fold5], m1_X[test_index_fold5], m1_y[train_index_fold5], m1_y[test_index_fold5]

############################## Method 2 ##############################

m2_X = m2.X
m2_y = m2.obs['disease__ontology_label']
kf.get_n_splits(m2_X)

## fold 1
xx = list(kf.split(m2_X))
xx0 = list(xx[0])

train_index_fold1 = xx0[0]
test_index_fold1 = xx0[1]

m2_X_train_fold1, m2_X_test_fold1, m2_y_train_fold1, m2_y_test_fold1 = m2_X[train_index_fold1], m2_X[test_index_fold1], m2_y[train_index_fold1], m2_y[test_index_fold1]

## fold 2
xx = list(kf.split(m2_X))
xx0 = list(xx[1])

train_index_fold2 = xx0[0]
test_index_fold2 = xx0[1]

m2_X_train_fold2, m2_X_test_fold2, m2_y_train_fold2, m2_y_test_fold2 = m2_X[train_index_fold2], m2_X[test_index_fold2], m2_y[train_index_fold2], m2_y[test_index_fold2]

## fold 3
xx = list(kf.split(m2_X))
xx0 = list(xx[2])

train_index_fold3 = xx0[0]
test_index_fold3 = xx0[1]

m2_X_train_fold3, m2_X_test_fold3, m2_y_train_fold3, m2_y_test_fold3 = m2_X[train_index_fold3], m2_X[test_index_fold3], m2_y[train_index_fold3], m2_y[test_index_fold3]

## fold 4
xx = list(kf.split(m2_X))
xx0 = list(xx[3])

train_index_fold4 = xx0[0]
test_index_fold4 = xx0[1]

m2_X_train_fold4, m2_X_test_fold4, m2_y_train_fold4, m2_y_test_fold4 = m2_X[train_index_fold4], m2_X[test_index_fold4], m2_y[train_index_fold4], m2_y[test_index_fold4]

## fold 5
xx = list(kf.split(m2_X))
xx0 = list(xx[4])

train_index_fold5 = xx0[0]
test_index_fold5 = xx0[1]

m2_X_train_fold5, m2_X_test_fold5, m2_y_train_fold5, m2_y_test_fold5 = m2_X[train_index_fold5], m2_X[test_index_fold5], m2_y[train_index_fold5], m2_y[test_index_fold5]

############################## Method 3 ##############################

m3_X = m3.X
m3_y = m3.obs['disease__ontology_label']
kf.get_n_splits(m3_X)

## fold 1
xx = list(kf.split(m3_X))
xx0 = list(xx[0])

train_index_fold1 = xx0[0]
test_index_fold1 = xx0[1]

m3_X_train_fold1, m3_X_test_fold1, m3_y_train_fold1, m3_y_test_fold1 = m3_X[train_index_fold1], m3_X[test_index_fold1], m3_y[train_index_fold1], m3_y[test_index_fold1]

## fold 2
xx = list(kf.split(m3_X))
xx0 = list(xx[1])

train_index_fold2 = xx0[0]
test_index_fold2 = xx0[1]

m3_X_train_fold2, m3_X_test_fold2, m3_y_train_fold2, m3_y_test_fold2 = m3_X[train_index_fold2], m3_X[test_index_fold2], m3_y[train_index_fold2], m3_y[test_index_fold2]

## fold 3
xx = list(kf.split(m3_X))
xx0 = list(xx[2])

train_index_fold3 = xx0[0]
test_index_fold3 = xx0[1]

m3_X_train_fold3, m3_X_test_fold3, m3_y_train_fold3, m3_y_test_fold3 = m3_X[train_index_fold3], m3_X[test_index_fold3], m3_y[train_index_fold3], m3_y[test_index_fold3]

## fold 4
xx = list(kf.split(m3_X))
xx0 = list(xx[3])

train_index_fold4 = xx0[0]
test_index_fold4 = xx0[1]

m3_X_train_fold4, m3_X_test_fold4, m3_y_train_fold4, m3_y_test_fold4 = m3_X[train_index_fold4], m3_X[test_index_fold4], m3_y[train_index_fold4], m3_y[test_index_fold4]

## fold 5
xx = list(kf.split(m3_X))
xx0 = list(xx[4])

train_index_fold5 = xx0[0]
test_index_fold5 = xx0[1]

m3_X_train_fold5, m3_X_test_fold5, m3_y_train_fold5, m3_y_test_fold5 = m3_X[train_index_fold5], m3_X[test_index_fold5], m3_y[train_index_fold5], m3_y[test_index_fold5]

############################## Method 4 ##############################

m4_X = m4.obsm["delta"]
m4_y = m4.obs['disease__ontology_label']
kf.get_n_splits(m4_X)

## fold 1
xx = list(kf.split(m4_X))
xx0 = list(xx[0])

train_index_fold1 = xx0[0]
test_index_fold1 = xx0[1]

m4_X_train_fold1, m4_X_test_fold1, m4_y_train_fold1, m4_y_test_fold1 = m4_X[train_index_fold1], m4_X[test_index_fold1], m4_y[train_index_fold1], m4_y[test_index_fold1]

## fold 2
xx = list(kf.split(m4_X))
xx0 = list(xx[1])

train_index_fold2 = xx0[0]
test_index_fold2 = xx0[1]

m4_X_train_fold2, m4_X_test_fold2, m4_y_train_fold2, m4_y_test_fold2 = m4_X[train_index_fold2], m4_X[test_index_fold2], m4_y[train_index_fold2], m4_y[test_index_fold2]

## fold 3
xx = list(kf.split(m4_X))
xx0 = list(xx[2])

train_index_fold3 = xx0[0]
test_index_fold3 = xx0[1]

m4_X_train_fold3, m4_X_test_fold3, m4_y_train_fold3, m4_y_test_fold3 = m4_X[train_index_fold3], m4_X[test_index_fold3], m4_y[train_index_fold3], m4_y[test_index_fold3]

## fold 4
xx = list(kf.split(m4_X))
xx0 = list(xx[3])

train_index_fold4 = xx0[0]
test_index_fold4 = xx0[1]

m4_X_train_fold4, m4_X_test_fold4, m4_y_train_fold4, m4_y_test_fold4 = m4_X[train_index_fold4], m4_X[test_index_fold4], m4_y[train_index_fold4], m4_y[test_index_fold4]

## fold 5
xx = list(kf.split(m4_X))
xx0 = list(xx[4])

train_index_fold5 = xx0[0]
test_index_fold5 = xx0[1]

m4_X_train_fold5, m4_X_test_fold5, m4_y_train_fold5, m4_y_test_fold5 = m4_X[train_index_fold5], m4_X[test_index_fold5], m4_y[train_index_fold5], m4_y[test_index_fold5]


###################################
# load test objects
###################################
## Load train data object
m1_test = sc.read_h5ad(par_m1_test)
m2_test = sc.read_h5ad(par_m2_test)
m3_test = sc.read_h5ad(par_m3_test)
m4_test = sc.read_h5ad(par_m4_test)

## fix disease ontology label for m4
m4_test.obs['disease__ontology_label'] = m4_test.obs['disease__ontology_label'].replace({'Parkinson disease': 1, 'normal': 0})

## define variables
X_test_m1 = m1_test.X
y_test_m1 = m1_test.obs['disease__ontology_label']
X_test_m2 = m2_test.X
y_test_m2 = m2_test.obs['disease__ontology_label']
X_test_m3 = m3_test.X
y_test_m3 = m3_test.obs['disease__ontology_label']
X_test_m4 = m4_test.obsm["delta"]
y_test_m4 = m4_test.obs['disease__ontology_label']

print(X_test_m1.shape, X_test_m2.shape, X_test_m3.shape, X_test_m4.shape)

num_genes = X_test_m1.shape[1]
num_PC = X_test_m2.shape[1]
num_NMF = X_test_m3.shape[1]
num_ETM = X_test_m4.shape[1]
print(num_genes, num_PC, num_NMF, num_ETM)

###################################
# Define functions
###################################

def dataset_split(X, y):
    # train-test split
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=5) #test is always the same
    # train-validation split
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.20)
    return X_train, X_val, y_train, y_val

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


##################################################################################################################################################################################################################
# Support vector machine
##################################################################################################################################################################################################################
###################################
# Produce SVM model 
###################################
## Method 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fold 1
SVM_m1_fold1 = svm.SVC(kernel='linear', C=1)
SVM_m1_fold1.fit(m1_X_train_fold1, m1_y_train_fold1)
SVM_m1_classifer_fold1 = SVM_m1_fold1.fit(m1_X_train_fold1, m1_y_train_fold1)
prediction_val, accuracy_val = validation(SVM_m1_classifer_fold1, m1_X_test_fold1, m1_y_test_fold1) 
## fold 2
SVM_m1_fold2 = svm.SVC(kernel='linear', C=1)
SVM_m1_fold2.fit(m1_X_train_fold2, m1_y_train_fold2)
SVM_m1_classifer_fold2 = SVM_m1_fold2.fit(m1_X_train_fold2, m1_y_train_fold2)
prediction_val, accuracy_val = validation(SVM_m1_classifer_fold2, m1_X_test_fold2, m1_y_test_fold2) 
## fold 3
SVM_m1_fold3 = svm.SVC(kernel='linear', C=1)
SVM_m1_fold3.fit(m1_X_train_fold3, m1_y_train_fold3)
SVM_m1_classifer_fold3 = SVM_m1_fold3.fit(m1_X_train_fold3, m1_y_train_fold3)
prediction_val, accuracy_val = validation(SVM_m1_classifer_fold3, m1_X_test_fold3, m1_y_test_fold3) 
## fold 4
SVM_m1_fold4 = svm.SVC(kernel='linear', C=1)
SVM_m1_fold4.fit(m1_X_train_fold4, m1_y_train_fold4)
SVM_m1_classifer_fold4 = SVM_m1_fold4.fit(m1_X_train_fold4, m1_y_train_fold4)
prediction_val, accuracy_val = validation(SVM_m1_classifer_fold4, m1_X_test_fold4, m1_y_test_fold4) 
## fold 5
SVM_m1_fold5 = svm.SVC(kernel='linear', C=1)
SVM_m1_fold5.fit(m1_X_train_fold5, m1_y_train_fold5)
SVM_m1_classifer_fold5 = SVM_m1_fold5.fit(m1_X_train_fold5, m1_y_train_fold5)
prediction_val, accuracy_val = validation(SVM_m1_classifer_fold5, m1_X_test_fold5, m1_y_test_fold5) 


## Method 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fold 1
SVM_m2_fold1 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m2_fold1.fit(m2_X_train_fold1, m2_y_train_fold1)
SVM_m2_classifer_fold1 = SVM_m2_fold1.fit(m2_X_train_fold1, m2_y_train_fold1)
prediction_val, accuracy_val = validation(SVM_m2_classifer_fold1, m2_X_test_fold1, m2_y_test_fold1) 
## fold 2
SVM_m2_fold2 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m2_fold2.fit(m2_X_train_fold2, m2_y_train_fold2)
SVM_m2_classifer_fold2 = SVM_m2_fold2.fit(m2_X_train_fold2, m2_y_train_fold2)
prediction_val, accuracy_val = validation(SVM_m2_classifer_fold2, m2_X_test_fold2, m2_y_test_fold2) 
## fold 3
SVM_m2_fold3 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m2_fold3.fit(m2_X_train_fold3, m2_y_train_fold3)
SVM_m2_classifer_fold3 = SVM_m2_fold3.fit(m2_X_train_fold3, m2_y_train_fold3)
prediction_val, accuracy_val = validation(SVM_m2_classifer_fold3, m2_X_test_fold3, m2_y_test_fold3) 
## fold 4
SVM_m2_fold4 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m2_fold4.fit(m2_X_train_fold4, m2_y_train_fold4)
SVM_m2_classifer_fold4 = SVM_m2_fold4.fit(m2_X_train_fold4, m2_y_train_fold4)
prediction_val, accuracy_val = validation(SVM_m2_classifer_fold4, m2_X_test_fold4, m2_y_test_fold4) 
## fold 5
SVM_m2_fold5 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m2_fold5.fit(m2_X_train_fold5, m2_y_train_fold5)
SVM_m2_classifer_fold5 = SVM_m2_fold5.fit(m2_X_train_fold5, m2_y_train_fold5)
prediction_val, accuracy_val = validation(SVM_m2_classifer_fold5, m2_X_test_fold5, m2_y_test_fold5) 


## Method 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fold 1
SVM_m3_fold1 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m3_fold1.fit(m3_X_train_fold1, m3_y_train_fold1)
SVM_m3_classifer_fold1 = SVM_m3_fold1.fit(m3_X_train_fold1, m3_y_train_fold1)
prediction_val, accuracy_val = validation(SVM_m3_classifer_fold1, m3_X_test_fold1, m3_y_test_fold1) 
## fold 2
SVM_m3_fold2 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m3_fold2.fit(m3_X_train_fold2, m3_y_train_fold2)
SVM_m3_classifer_fold2 = SVM_m3_fold2.fit(m3_X_train_fold2, m3_y_train_fold2)
prediction_val, accuracy_val = validation(SVM_m3_classifer_fold2, m3_X_test_fold2, m3_y_test_fold2) 
## fold 3
SVM_m3_fold3 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m3_fold3.fit(m3_X_train_fold3, m3_y_train_fold3)
SVM_m3_classifer_fold3 = SVM_m3_fold3.fit(m3_X_train_fold3, m3_y_train_fold3)
prediction_val, accuracy_val = validation(SVM_m3_classifer_fold3, m3_X_test_fold3, m3_y_test_fold3) 
## fold 4
SVM_m3_fold4 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m3_fold4.fit(m3_X_train_fold4, m3_y_train_fold4)
SVM_m3_classifer_fold4 = SVM_m3_fold4.fit(m3_X_train_fold4, m3_y_train_fold4)
prediction_val, accuracy_val = validation(SVM_m3_classifer_fold4, m3_X_test_fold4, m3_y_test_fold4) 
## fold 5
SVM_m3_fold5 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m3_fold5.fit(m3_X_train_fold5, m3_y_train_fold5)
SVM_m3_classifer_fold5 = SVM_m3_fold5.fit(m3_X_train_fold5, m3_y_train_fold5)
prediction_val, accuracy_val = validation(SVM_m3_classifer_fold5, m3_X_test_fold5, m3_y_test_fold5) 

## Method 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fold 1
SVM_m4_fold1 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m4_fold1.fit(m4_X_train_fold1, m4_y_train_fold1)
SVM_m4_classifer_fold1 = SVM_m4_fold1.fit(m4_X_train_fold1, m4_y_train_fold1)
prediction_val, accuracy_val = validation(SVM_m4_classifer_fold1, m4_X_test_fold1, m4_y_test_fold1) 
## fold 2
SVM_m4_fold2 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m4_fold2.fit(m4_X_train_fold2, m4_y_train_fold2)
SVM_m4_classifer_fold2 = SVM_m4_fold2.fit(m4_X_train_fold2, m4_y_train_fold2)
prediction_val, accuracy_val = validation(SVM_m4_classifer_fold2, m4_X_test_fold2, m4_y_test_fold2) 
## fold 3
SVM_m4_fold3 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m4_fold3.fit(m4_X_train_fold3, m4_y_train_fold3)
SVM_m4_classifer_fold3 = SVM_m4_fold3.fit(m4_X_train_fold3, m4_y_train_fold3)
prediction_val, accuracy_val = validation(SVM_m4_classifer_fold3, m4_X_test_fold3, m4_y_test_fold3) 
## fold 4
SVM_m4_fold4 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m4_fold4.fit(m4_X_train_fold4, m4_y_train_fold4)
SVM_m4_classifer_fold4 = SVM_m4_fold4.fit(m4_X_train_fold4, m4_y_train_fold4)
prediction_val, accuracy_val = validation(SVM_m4_classifer_fold4, m4_X_test_fold4, m4_y_test_fold4) 
## fold 5
SVM_m4_fold5 = svm.SVC(kernel='linear', C=1, probability=True)
SVM_m4_fold5.fit(m4_X_train_fold5, m4_y_train_fold5)
SVM_m4_classifer_fold5 = SVM_m4_fold5.fit(m4_X_train_fold5, m4_y_train_fold5)
prediction_val, accuracy_val = validation(SVM_m4_classifer_fold5, m4_X_test_fold5, m4_y_test_fold5) 


###################################
# Test 
###################################

############################## Method 1 ##############################
##### fold 1 #####
## prediction
prediction_test = SVM_m1_fold1.predict(X_test_m1)
#predictions = SVM_m1_classifer_fold1.predict_proba(X_test_m1)
## accuracy
accuracy_test_m1_fold1 = metrics.accuracy_score(y_test_m1, prediction_test)
## BA
BA_test_m1_fold1 = balanced_accuracy_score(y_test_m1, prediction_test)
## F1
F1_test_m1_fold1 = f1_score(y_test_m1, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m1, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m1','m1','m1','m1']
}
confusion_matrix_m1_fold1 = pd.DataFrame(data) 

##### fold 2 #####
## prediction
prediction_test = SVM_m1_fold2.predict(X_test_m1)
#predictions = SVM_m1_classifer_fold2.predict_proba(X_test_m1)
## accuracy
accuracy_test_m1_fold2 = metrics.accuracy_score(y_test_m1, prediction_test)
## BA
BA_test_m1_fold2 = balanced_accuracy_score(y_test_m1, prediction_test)
## F1
F1_test_m1_fold2 = f1_score(y_test_m1, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m1, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m1','m1','m1','m1']
}
confusion_matrix_m1_fold2 = pd.DataFrame(data) 

##### fold 3 #####
## prediction
prediction_test = SVM_m1_fold3.predict(X_test_m1)
#predictions = SVM_m1_classifer_fold3.predict_proba(X_test_m1)
## accuracy
accuracy_test_m1_fold3 = metrics.accuracy_score(y_test_m1, prediction_test)
## BA
BA_test_m1_fold3 = balanced_accuracy_score(y_test_m1, prediction_test)
## F1
F1_test_m1_fold3 = f1_score(y_test_m1, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m1, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m1','m1','m1','m1']
}
confusion_matrix_m1_fold3 = pd.DataFrame(data) 

##### fold 4 #####
## prediction
prediction_test = SVM_m1_fold4.predict(X_test_m1)
#predictions = SVM_m1_classifer_fold4.predict_proba(X_test_m1)
## accuracy
accuracy_test_m1_fold4 = metrics.accuracy_score(y_test_m1, prediction_test)
## BA
BA_test_m1_fold4 = balanced_accuracy_score(y_test_m1, prediction_test)
## F1
F1_test_m1_fold4 = f1_score(y_test_m1, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m1, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m1','m1','m1','m1']
}
confusion_matrix_m1_fold4 = pd.DataFrame(data) 

##### fold 5 #####
## prediction
prediction_test = SVM_m1_fold5.predict(X_test_m1)
#predictions = SVM_m1_classifer_fold5.predict_proba(X_test_m1)
## accuracy
accuracy_test_m1_fold5 = metrics.accuracy_score(y_test_m1, prediction_test)
## BA
BA_test_m1_fold5 = balanced_accuracy_score(y_test_m1, prediction_test)
## F1
F1_test_m1_fold5 = f1_score(y_test_m1, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m1, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m1','m1','m1','m1']
}
confusion_matrix_m1_fold5 = pd.DataFrame(data) 


############################## Method 2 ##############################
##### fold 1 #####
## prediction
prediction_test = SVM_m2_fold1.predict(X_test_m2)
predictions = SVM_m2_classifer_fold1.predict_proba(X_test_m2)
## accuracy
accuracy_test_m2_fold1 = metrics.accuracy_score(y_test_m2, prediction_test)
## BA
BA_test_m2_fold1 = balanced_accuracy_score(y_test_m2, prediction_test)
## F1
F1_test_m2_fold1 = f1_score(y_test_m2, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m2, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m2','m2','m2','m2']
}
confusion_matrix_m2_fold1 = pd.DataFrame(data) 


##### fold 2 #####
## prediction
prediction_test = SVM_m2_fold2.predict(X_test_m2)
predictions = SVM_m2_classifer_fold2.predict_proba(X_test_m2)
## accuracy
accuracy_test_m2_fold2 = metrics.accuracy_score(y_test_m2, prediction_test)
## BA
BA_test_m2_fold2 = balanced_accuracy_score(y_test_m2, prediction_test)
## F1
F1_test_m2_fold2 = f1_score(y_test_m2, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m2, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m2','m2','m2','m2']
}
confusion_matrix_m2_fold2 = pd.DataFrame(data) 


##### fold 3 #####
## prediction
prediction_test = SVM_m2_fold3.predict(X_test_m2)
predictions = SVM_m2_classifer_fold3.predict_proba(X_test_m2)
## accuracy
accuracy_test_m2_fold3 = metrics.accuracy_score(y_test_m2, prediction_test)
## BA
BA_test_m2_fold3 = balanced_accuracy_score(y_test_m2, prediction_test)
## F1
F1_test_m2_fold3 = f1_score(y_test_m2, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m2, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m2','m2','m2','m2']
}
confusion_matrix_m2_fold3 = pd.DataFrame(data) 


##### fold 4 #####
## prediction
prediction_test = SVM_m2_fold4.predict(X_test_m2)
predictions = SVM_m2_classifer_fold4.predict_proba(X_test_m2)
## accuracy
accuracy_test_m2_fold4 = metrics.accuracy_score(y_test_m2, prediction_test)
## BA
BA_test_m2_fold4 = balanced_accuracy_score(y_test_m2, prediction_test)
## F1
F1_test_m2_fold4 = f1_score(y_test_m2, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m2, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m2','m2','m2','m2']
}
confusion_matrix_m2_fold4 = pd.DataFrame(data) 

##### fold 5 #####
## prediction
prediction_test = SVM_m2_fold5.predict(X_test_m2)
predictions = SVM_m2_classifer_fold5.predict_proba(X_test_m2)
## accuracy
accuracy_test_m2_fold5 = metrics.accuracy_score(y_test_m2, prediction_test)
## BA
BA_test_m2_fold5 = balanced_accuracy_score(y_test_m2, prediction_test)
## F1
F1_test_m2_fold5 = f1_score(y_test_m2, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m2, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m2','m2','m2','m2']
}
confusion_matrix_m2_fold5 = pd.DataFrame(data) 


############################## Method 3 ##############################
##### fold 1 #####
## prediction
prediction_test = SVM_m3_fold1.predict(X_test_m3)
predictions = SVM_m3_classifer_fold1.predict_proba(X_test_m3)

## accuracy
accuracy_test_m3_fold1 = metrics.accuracy_score(y_test_m3, prediction_test)
## BA
BA_test_m3_fold1 = balanced_accuracy_score(y_test_m3, prediction_test)
## F1
F1_test_m3_fold1 = f1_score(y_test_m3, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m3, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m3','m3','m3','m3']
}
confusion_matrix_m3_fold1 = pd.DataFrame(data) 

#need to get the index
df_index = pd.DataFrame(y_test_m3)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(SVM_par_cell_index_m3_fold1)

##### fold 2 #####
## prediction
prediction_test = SVM_m3_fold2.predict(X_test_m3)
predictions = SVM_m3_classifer_fold2.predict_proba(X_test_m3)

## accuracy
accuracy_test_m3_fold2 = metrics.accuracy_score(y_test_m3, prediction_test)
## BA
BA_test_m3_fold2 = balanced_accuracy_score(y_test_m3, prediction_test)
## F1
F1_test_m3_fold2 = f1_score(y_test_m3, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m3, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m3','m3','m3','m3']
}
confusion_matrix_m3_fold2 = pd.DataFrame(data) 

#need to get the index
df_index = pd.DataFrame(y_test_m3)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(SVM_par_cell_index_m3_fold2)

##### fold 3 #####
## prediction
prediction_test = SVM_m3_fold3.predict(X_test_m3)
predictions = SVM_m3_classifer_fold3.predict_proba(X_test_m3)

## accuracy
accuracy_test_m3_fold3 = metrics.accuracy_score(y_test_m3, prediction_test)
## BA
BA_test_m3_fold3 = balanced_accuracy_score(y_test_m3, prediction_test)
## F1
F1_test_m3_fold3 = f1_score(y_test_m3, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m3, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m3','m3','m3','m3']
}
confusion_matrix_m3_fold3 = pd.DataFrame(data) 

##### fold 4 #####
## prediction
prediction_test = SVM_m3_fold4.predict(X_test_m3)
predictions = SVM_m3_classifer_fold4.predict_proba(X_test_m3)

## accuracy
accuracy_test_m3_fold4 = metrics.accuracy_score(y_test_m3, prediction_test)
## BA
BA_test_m3_fold4 = balanced_accuracy_score(y_test_m3, prediction_test)
## F1
F1_test_m3_fold4 = f1_score(y_test_m3, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m3, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m3','m3','m3','m3']
}
confusion_matrix_m3_fold4 = pd.DataFrame(data) 

##### fold 5 #####
## prediction
prediction_test = SVM_m3_fold5.predict(X_test_m3)
predictions = SVM_m3_classifer_fold5.predict_proba(X_test_m3)

## accuracy
accuracy_test_m3_fold5 = metrics.accuracy_score(y_test_m3, prediction_test)
## BA
BA_test_m3_fold5 = balanced_accuracy_score(y_test_m3, prediction_test)
## F1
F1_test_m3_fold5 = f1_score(y_test_m3, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m3, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m3','m3','m3','m3']
}
confusion_matrix_m3_fold5 = pd.DataFrame(data) 


############################## Method 4 ##############################
##### fold 1 #####
## prediction
prediction_test = SVM_m4_fold1.predict(X_test_m4)
predictions = SVM_m4_classifer_fold1.predict_proba(X_test_m4)

## accuracy
accuracy_test_m4_fold1 = metrics.accuracy_score(y_test_m4, prediction_test)
## BA
BA_test_m4_fold1 = balanced_accuracy_score(y_test_m4, prediction_test)
## F1
F1_test_m4_fold1 = f1_score(y_test_m4, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m4, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m4','m4','m4','m4']
}
confusion_matrix_m4_fold1 = pd.DataFrame(data) 


##### fold 2 #####
## prediction
prediction_test = SVM_m4_fold2.predict(X_test_m4)
predictions = SVM_m4_classifer_fold2.predict_proba(X_test_m4)

## accuracy
accuracy_test_m4_fold2 = metrics.accuracy_score(y_test_m4, prediction_test)
## BA
BA_test_m4_fold2 = balanced_accuracy_score(y_test_m4, prediction_test)
## F1
F1_test_m4_fold2 = f1_score(y_test_m4, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m4, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m4','m4','m4','m4']
}
confusion_matrix_m4_fold2 = pd.DataFrame(data) 


##### fold 3 #####
## prediction
prediction_test = SVM_m4_fold3.predict(X_test_m4)
predictions = SVM_m4_classifer_fold3.predict_proba(X_test_m4)

## accuracy
accuracy_test_m4_fold3 = metrics.accuracy_score(y_test_m4, prediction_test)
## BA
BA_test_m4_fold3 = balanced_accuracy_score(y_test_m4, prediction_test)
## F1
F1_test_m4_fold3 = f1_score(y_test_m4, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m4, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m4','m4','m4','m4']
}
confusion_matrix_m4_fold3 = pd.DataFrame(data) 

##### fold 4 #####
## prediction
prediction_test = SVM_m4_fold4.predict(X_test_m4)
predictions = SVM_m4_classifer_fold4.predict_proba(X_test_m4)

## accuracy
accuracy_test_m4_fold4 = metrics.accuracy_score(y_test_m4, prediction_test)
## BA
BA_test_m4_fold4 = balanced_accuracy_score(y_test_m4, prediction_test)
## F1
F1_test_m4_fold4 = f1_score(y_test_m4, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m4, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m4','m4','m4','m4']
}
confusion_matrix_m4_fold4 = pd.DataFrame(data) 


##### fold 5 #####
## prediction
prediction_test = SVM_m4_fold5.predict(X_test_m4)
predictions = SVM_m4_classifer_fold5.predict_proba(X_test_m4)

## accuracy
accuracy_test_m4_fold5 = metrics.accuracy_score(y_test_m4, prediction_test)
## BA
BA_test_m4_fold5 = balanced_accuracy_score(y_test_m4, prediction_test)
## F1
F1_test_m4_fold5 = f1_score(y_test_m4, prediction_test, pos_label=1)
# confusuon matrix values
cm_array = confusion_matrix(y_test_m4, prediction_test)
# save the confusion matrix
data = {'counts': [cm_array[1,1], cm_array[0,0], cm_array[1,0], cm_array[0,1]],
'true_label': ['PD', 'Healthy','PD', 'Healthy'],
'predict_label': ['PD', 'Healthy', 'Healthy', 'PD' ],
'method': ['m4','m4','m4','m4']
}
confusion_matrix_m4_fold5 = pd.DataFrame(data) 

## Combine all results for test run
test_results = [accuracy_test_m1_fold1, accuracy_test_m1_fold2, accuracy_test_m1_fold3, accuracy_test_m1_fold4, accuracy_test_m1_fold5,
accuracy_test_m2_fold1,accuracy_test_m2_fold2,accuracy_test_m2_fold3,accuracy_test_m2_fold4,accuracy_test_m2_fold5, 
accuracy_test_m3_fold1,accuracy_test_m3_fold2,accuracy_test_m3_fold3,accuracy_test_m3_fold4,accuracy_test_m3_fold5,
accuracy_test_m4_fold1,accuracy_test_m4_fold2,accuracy_test_m4_fold3,accuracy_test_m4_fold4,accuracy_test_m4_fold5,
BA_test_m1_fold1, BA_test_m1_fold2, BA_test_m1_fold3, BA_test_m1_fold4, BA_test_m1_fold5,
BA_test_m2_fold1,BA_test_m2_fold2,BA_test_m2_fold3,BA_test_m2_fold4,BA_test_m2_fold5, 
BA_test_m3_fold1,BA_test_m3_fold2,BA_test_m3_fold3,BA_test_m3_fold4,BA_test_m3_fold5,
BA_test_m4_fold1,BA_test_m4_fold2,BA_test_m4_fold3,BA_test_m4_fold4,BA_test_m4_fold5,
F1_test_m1_fold1, F1_test_m1_fold2, F1_test_m1_fold3, F1_test_m1_fold4, F1_test_m1_fold5,
F1_test_m2_fold1,F1_test_m2_fold2,F1_test_m2_fold3,F1_test_m2_fold4,F1_test_m2_fold5, 
F1_test_m3_fold1,F1_test_m3_fold2,F1_test_m3_fold3,F1_test_m3_fold4,F1_test_m3_fold5,
F1_test_m4_fold1,F1_test_m4_fold2,F1_test_m4_fold3,F1_test_m4_fold4,F1_test_m4_fold5]

df = pd.DataFrame(test_results)
df
df2 = df.assign( metric = ['accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy','accuracy',
'BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA','BA',
'F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1','F1'])
df2 = df2.assign( method = ['m1','m1','m1','m1','m1',
'm2','m2','m2','m2','m2',
'm3','m3','m3','m3','m3',
'm4','m4','m4','m4','m4',
'm1','m1','m1','m1','m1',
'm2','m2','m2','m2','m2',
'm3','m3','m3','m3','m3',
'm4','m4','m4','m4','m4',
'm1','m1','m1','m1','m1',
'm2','m2','m2','m2','m2',
'm3','m3','m3','m3','m3',
'm4','m4','m4','m4','m4'])

df2.to_csv(par_metric_summary_frame_SVM)

## create Confusion matrix output dataframe
df = pd.concat([confusion_matrix_m1_fold1,confusion_matrix_m1_fold2,confusion_matrix_m1_fold3,confusion_matrix_m1_fold4,confusion_matrix_m1_fold5,
confusion_matrix_m2_fold1,confusion_matrix_m2_fold2,confusion_matrix_m2_fold3,confusion_matrix_m2_fold4,confusion_matrix_m2_fold5,
confusion_matrix_m3_fold1,confusion_matrix_m3_fold2,confusion_matrix_m3_fold3,confusion_matrix_m3_fold4,confusion_matrix_m3_fold5,
confusion_matrix_m4_fold1,confusion_matrix_m4_fold2,confusion_matrix_m4_fold3,confusion_matrix_m4_fold4,confusion_matrix_m4_fold5])

df

df.to_csv(par_confusion_matrix_combo_SVM)


