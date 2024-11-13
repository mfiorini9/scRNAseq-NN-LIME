# scRNAseq-NN-LIME

This repository contains the code to train a neural network for classifying individual cells based on their treatment or disease status, using single-cell or single-nucleus RNA sequencing data (scnRNAseq). Additionally, it includes the code to decode the trained models using Local Interpretable Model-Agnostic Explanations (LIME) to identify the most important genes that distinguish treated or diseased cells.

We have demonstrated the application of this framework in our manuscript titled "Interpretable Machine Learning Classifiers Implicate GPC6 in Parkinson’s Disease from Single Nuclei Midbrain Transcriptomes." The preprint is available on bioRxiv.

## Contents
-  [Introduction](#introduction)
-  [How to run](#how-to-run)

---

## Introduction
Supervised machine learning (ML) models can be implemented to classify samples based on treatment of disease status. In the context of scnRNAseq, ML classifiers can be trained to accurately distinguish between the transcriptomes of treated/diseased and untreated/healthy cells, which would open avenues for deeper investigation into the underlying “reasons” for a particular cell receiving a particular classification. However, many ML classifiers operate as “black boxes”, making it difficult to interpret their decision-making processes. To address this, interpretable explainers like LIME can be applied to reveal the most influential features driving a classification decision. For scnRNAseq data, these influential features would constitute the expression profiles of specific genes. We have demonstrated that LIME-identified genes are not only valuable for in silico classifiation but also represent relevant biological perturbations. Furthermore, we have demonstrated the pan-dataset generalizability of LIME-identified genes, reinforcing that they capture replicable readouts across seuencign experiments. For more information please see our pre-print.

## How to run
We have assembled a python script for users that wish to apply this framework to their own analyses. To access this code, please see the [NN-LIME](https://github.com/mfiorini9/scRNAseq-NN-LIME/tree/main/NN-LIME) folder in this repository.

### License
This project is licensed under the MIT License. See the [LICENSE](https://github.com/mfiorini9/scRNAseq-NN-LIME/blob/main/LICENSE) file for details.

## Acknowledgement
scRNAseq-NN-LIME was produced by Michael Fiorini with associate of Jialun Li, Edward Fon, Sali Farhan, and Rhalena Thomas. 
