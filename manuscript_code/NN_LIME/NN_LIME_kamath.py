## This script was used to test/train the NN models and decode them with LIME for the Kamath object. We use astrocytes as an example. 

salloc -A def-sfarhan --time=0-5 -c 1 --mem=300g

module load StdEnv/2020 

module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

###################################
# Adjustable parameters
###################################
## all methods
par_m1_test = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/HVG_test_kamath_astro.h5ad'
par_m1_train = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/HVG_train_kamath_astro.h5ad'

DNN_par_cell_index_m1_fold1 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold1_kamath_astro.csv'
DNN_par_LIME_out_m1_fold1 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold1_kamath_astro.csv'

DNN_par_cell_index_m1_fold2 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold2_kamath_astro.csv'
DNN_par_LIME_out_m1_fold2 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold2_kamath_astro.csv'

DNN_par_cell_index_m1_fold3 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold3_kamath_astro.csv'
DNN_par_LIME_out_m1_fold3 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold3_kamath_astro.csv'

DNN_par_cell_index_m1_fold4 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold4_kamath_astro.csv'
DNN_par_LIME_out_m1_fold4 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold4_kamath_astro.csv'

DNN_par_cell_index_m1_fold5 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold5_kamath_astro.csv'
DNN_par_LIME_out_m1_fold5 = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold5_kamath_astro.csv'

par_metric_summary_frame_DNN = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_performance_summary_kamath_astro.csv'
par_confusion_matrix_combo_DNN = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_CM_kamath_astro.csv'


###################################
# load library
###################################
## DNN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn import metrics
import sklearn.metrics as metrics
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn import svm
import copy
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, f1_score, balanced_accuracy_score, make_scorer
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from math import sqrt
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report,confusion_matrix
import random
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, f1_score, balanced_accuracy_score
import lime
import lime.lime_tabular

from pprint import pprint
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.metrics import f1_score

from sklearn.metrics import make_scorer, f1_score

from sklearn.model_selection import cross_val_score, KFold
from sklearn.inspection import permutation_importance
import time

import lime
import lime.lime_tabular

from sklearn import metrics
from sklearn.metrics import roc_auc_score

from sklearn.model_selection import KFold

###################################
# preprocessing adata objects
###################################
## Load train data object
m1 = sc.read_h5ad(par_m1_train)

###################################
# 5-fold cross validation for each method
###################################
kf = KFold(n_splits=5)

## set variables
m1_X = m1.X
m1_y = m1.obs['Disease_Status']
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


###################################
# load test objects
###################################
## Load train data object
m1_test = sc.read_h5ad(par_m1_test)

## define variables
X_test_m1 = m1_test.X
y_test_m1 = m1_test.obs['Disease_Status']

print(X_test_m1.shape)

num_genes = X_test_m1.shape[1]
print(num_genes)

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

###################################
# Produce DNN model 
###################################
## Method 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fold1
MLP_m1_fold1 = MLPClassifier( activation='relu', solver='adam', max_iter=500)
MLP_m1_fold1.fit(m1_X_train_fold1, m1_y_train_fold1)
MLP_m1_classifier_fold1 = MLP_m1_fold1.fit(m1_X_train_fold1, m1_y_train_fold1)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold1, m1_X_test_fold1, m1_y_test_fold1) 
## fold2
MLP_m1_fold2 = MLPClassifier( activation='relu', solver='adam', max_iter=500)
MLP_m1_fold2.fit(m1_X_train_fold2, m1_y_train_fold2)
MLP_m1_classifier_fold2 = MLP_m1_fold2.fit(m1_X_train_fold2, m1_y_train_fold2)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold2, m1_X_test_fold2, m1_y_test_fold2) 
## fold3
MLP_m1_fold3 = MLPClassifier( activation='relu', solver='adam', max_iter=500)
MLP_m1_fold3.fit(m1_X_train_fold3, m1_y_train_fold3)
MLP_m1_classifier_fold3 = MLP_m1_fold3.fit(m1_X_train_fold3, m1_y_train_fold3)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold3, m1_X_test_fold3, m1_y_test_fold3) 
## fold4
MLP_m1_fold4 = MLPClassifier( activation='relu', solver='adam', max_iter=500)
MLP_m1_fold4.fit(m1_X_train_fold4, m1_y_train_fold4)
MLP_m1_classifier_fold4 = MLP_m1_fold4.fit(m1_X_train_fold4, m1_y_train_fold4)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold4, m1_X_test_fold4, m1_y_test_fold4) 
## fold5
MLP_m1_fold5 = MLPClassifier( activation='relu', solver='adam', max_iter=500)
MLP_m1_fold5.fit(m1_X_train_fold5, m1_y_train_fold5)
MLP_m1_classifier_fold5 = MLP_m1_fold5.fit(m1_X_train_fold5, m1_y_train_fold5)
prediction_val, accuracy_val = validation(MLP_m1_classifier_fold5, m1_X_test_fold5, m1_y_test_fold5) 

###################################
# Test
## NOTE: Add F1 score and balanced accuracy to the predict function
## NOTE: Add code here to get the confusion matrix values 
###################################

############################## Method 1 ##############################
##### fold1 #####
## prediction
prediction_test = MLP_m1_fold1.predict(X_test_m1)
predictions = MLP_m1_classifier_fold1.predict_proba(X_test_m1)

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

#need to get the index
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(DNN_par_cell_index_m1_fold1)


####### LIME #######
## convert array to dataframe
df=pd.DataFrame(m1_X_train_fold1,) 
df_test=pd.DataFrame(X_test_m1,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)


explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['Normal', 'Parkinson disease'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

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

## write LIME output csv
df.to_csv(DNN_par_LIME_out_m1_fold1)


##### fold2 #####
## prediction
prediction_test = MLP_m1_fold2.predict(X_test_m1)
predictions = MLP_m1_classifier_fold2.predict_proba(X_test_m1)

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

#need to get the index
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(DNN_par_cell_index_m1_fold2)


####### LIME #######
## convert array to dataframe
df=pd.DataFrame(m1_X_train_fold2,) 
df_test=pd.DataFrame(X_test_m1,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## LIME has one explainer for all the models
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['Normal', 'Parkinson disease'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## modify to concatonate to an empty list.
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

## write LIME output csv
df.to_csv(DNN_par_LIME_out_m1_fold2)


##### fold3 #####
## prediction
prediction_test = MLP_m1_fold3.predict(X_test_m1)
predictions = MLP_m1_classifier_fold3.predict_proba(X_test_m1)

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

#need to get the index
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(DNN_par_cell_index_m1_fold3)


####### LIME #######
## convert array to dataframe
df=pd.DataFrame(m1_X_train_fold3,) 
df_test=pd.DataFrame(X_test_m1,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## LIME has one explainer for all the models
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['Normal', 'Parkinson disease'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## modify to concatonate to an empty list.
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

## write LIME output csv
df.to_csv(DNN_par_LIME_out_m1_fold3)


##### fold4 #####
## prediction
prediction_test = MLP_m1_fold4.predict(X_test_m1)
predictions = MLP_m1_classifier_fold4.predict_proba(X_test_m1)

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

#need to get the index
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(DNN_par_cell_index_m1_fold4)


####### LIME #######
## convert array to dataframe
df=pd.DataFrame(m1_X_train_fold4,) 
df_test=pd.DataFrame(X_test_m1,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## LIME has one explainer for all the models
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['Normal', 'Parkinson disease'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## modify to concatonate to an empty list.
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

## write LIME output csv
df.to_csv(DNN_par_LIME_out_m1_fold4)



##### fold5 #####
## prediction
prediction_test = MLP_m1_fold5.predict(X_test_m1)
predictions = MLP_m1_classifier_fold5.predict_proba(X_test_m1)

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

#need to get the index
df_index = pd.DataFrame(y_test_m1)
df_index['predicted_label'] = prediction_test
df_index['cell_index'] = df_index.index

# export to csv file
df_index.to_csv(DNN_par_cell_index_m1_fold5)


####### LIME #######
## convert array to dataframe
df=pd.DataFrame(m1_X_train_fold5,) 
df_test=pd.DataFrame(X_test_m1,)
df_test.index = df_index['cell_index']
num_cells = np.shape(df_test.values)[0]  
print(num_cells)

## LIME has one explainer for all the models
explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['Normal', 'Parkinson disease'], feature_names=m1_test.var_names.tolist(), verbose=True, mode='classification')

## modify to concatonate to an empty list.
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

## write LIME output csv
df.to_csv(DNN_par_LIME_out_m1_fold5)

## Combine all results for test run
test_results = [accuracy_test_m1_fold1, accuracy_test_m1_fold2, accuracy_test_m1_fold3, accuracy_test_m1_fold4, accuracy_test_m1_fold5,
                BA_test_m1_fold1, BA_test_m1_fold2, BA_test_m1_fold3, BA_test_m1_fold4, BA_test_m1_fold5,
                F1_test_m1_fold1, F1_test_m1_fold2, F1_test_m1_fold3, F1_test_m1_fold4, F1_test_m1_fold5]

df = pd.DataFrame(test_results)
df
df2 = df.assign( metric = ['accuracy','accuracy','accuracy','accuracy','accuracy',
                            'BA','BA','BA','BA','BA',
                            'F1','F1','F1','F1','F1'])
df2 = df2.assign( method = ['m1','m1','m1','m1','m1',
                            'm1','m1','m1','m1','m1',
                            'm1','m1','m1','m1','m1'])

df2.to_csv(par_metric_summary_frame_DNN)

## create Confusion matrix output dataframe
df = pd.concat([confusion_matrix_m1_fold1,confusion_matrix_m1_fold2,confusion_matrix_m1_fold3,confusion_matrix_m1_fold4,confusion_matrix_m1_fold5])
df

df.to_csv(par_confusion_matrix_combo_DNN)



