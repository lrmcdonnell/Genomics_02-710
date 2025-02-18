# classification.py
# HW2, Computational Genomics, Spring 2025
# andrewid: lmcdonne

# WARNING: Do not change the file name; Autograder expects it.

import sys

import numpy as np
from scipy.sparse import csc_matrix, save_npz, load_npz

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestClassifier

import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

def get_top_gene_filter(data, n_keep = 2000):
    """Select top n_keep most dispersed genes.

    Args:
        data (n x m matrix): input gene expression data of shape num_cells x num_genes
        n_keep (int): number of genes to be kepted after filtration; default 2000

    Returns:
        filter (array of length n_keep): an array of column indices that can be used as an
            index to keep only certain genes in data. Each element of filter is the column
            index of a highly-dispersed gene in data.
    """
    means=np.mean(data, axis=0) 
    means=np.where(means != 0, means, 0.000001)
    vars=np.var(data, axis=0)
    disps=vars/means
    idx=disps.argsort()[-n_keep:] 
    return idx

def reduce_dimensionality_pca(filtered_train_gene_expression, filtered_test_gene_expression, n_components = 20):
    """Train a PCA model and use it to reduce the training and testing data.
    
    Args:
        filtered_train_gene_expression (n_train x num_top_genes matrix): input filtered training expression data 
        filtered_test_gene_expression (n_test x num_top_genes matrix): input filtered test expression data 
        
    Return:
        (reduced_train_data, reduced_test_data): a tuple of
            1. The filtered training data transformed to the PC space.
            2. The filtered test data transformed to the PC space.
    """
    combined_data = np.concatenate([filtered_train_gene_expression, filtered_test_gene_expression], axis=0)
    
    pca = PCA(n_components=n_components)
    pca.fit(combined_data)
    
    reduced_train_data = pca.transform(filtered_train_gene_expression)
    reduced_test_data = pca.transform(filtered_test_gene_expression)
    
    return reduced_train_data, reduced_test_data

def plot_transformed_cells(reduced_train_data, train_labels):
    """Plot the PCA-reduced training data using just the first 2 principal components.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        None

    """
    df = pd.DataFrame(reduced_train_data, columns=[f'PC{i+1}' for i in range(reduced_train_data.shape[1])])
    df['cell_type'] = train_labels

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', hue='cell_type', data=df, palette='tab10', alpha=0.7)
    plt.title("PCA Projection of Training Cells")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(title='Cell Type')
    plt.tight_layout()
    plt.show()
    
def train_and_evaluate_rf_classifier(reduced_train_data, reduced_test_data, train_labels, test_labels):
    """Train and evaluate a simple Random Forest classification pipeline.
    
    Before passing the data to the RF module, this function scales the data such that the mean
    is 0 and the variance is 1.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        (classifier, score): a tuple consisting of
            1. classifier: the trained classifier
            2. The score (accuracy) of the classifier on the test data.

    """
    scaler = StandardScaler()
    scaled_train = scaler.fit_transform(reduced_train_data)
    scaled_test = scaler.transform(reduced_test_data)
    rf_classifier = RandomForestClassifier(random_state=42)
    rf_classifier.fit(scaled_train, train_labels)
    
    score = rf_classifier.score(scaled_test, test_labels)
    print("Random Forest Test Accuracy: {:.2f}%".format(score * 100))
    
    return rf_classifier, score
        
if __name__ == "__main__":
    train_gene_expression = np.load(sys.argv[1])['train']
    test_gene_expression = np.load(sys.argv[2])['test']
    train_labels = np.load(sys.argv[3])
    test_labels = np.load(sys.argv[4])
    
    top_gene_filter = get_top_gene_filter(train_gene_expression)
    filtered_test_gene_expression = test_gene_expression[:, top_gene_filter]
    filtered_train_gene_expression = train_gene_expression[:, top_gene_filter]
        
    mode = sys.argv[5]
    if mode == "rf_pipeline":
        (reduced_train_data,reduced_test_data) = reduce_dimensionality_pca(filtered_train_gene_expression, filtered_test_gene_expression, n_components = 20)
        plot_transformed_cells(reduced_train_data, train_labels)
        train_and_evaluate_rf_classifier(reduced_train_data, reduced_test_data, train_labels, test_labels)