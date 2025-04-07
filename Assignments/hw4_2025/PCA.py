import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path):
    try:
        data = pd.read_csv(file_path)
        print(f"Data loaded successfully with shape: {data.shape}")
        
        X = data.iloc[:, :-1]
        y = data.iloc[:, -1]
        
        print(f"Number of samples: {X.shape[0]}")
        print(f"Number of genes: {X.shape[1]}")
        print(f"Classes: {y.unique()}")
        
        return X, y
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None

def standardize_data(data):
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(data)
    return standardized_data, scaler

def apply_pca(standardized_data, n_components=2):
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(standardized_data)
    explained_variance_ratio = pca.explained_variance_ratio_
    return principal_components, pca, explained_variance_ratio

def plot_pca_results(principal_components, y, explained_variance_ratio):
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
    pca_df['Class'] = y
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Class', style='Class')
    plt.title('PCA of Gene Expression Data')
    plt.xlabel(f'PC1 ({explained_variance_ratio[0]*100:.2f}%)')
    plt.ylabel(f'PC2 ({explained_variance_ratio[1]*100:.2f}%)')
    
    plt.savefig('figures/pca_plot.png')
    plt.close()

def main():
    file_path = "provided_data/gene_expression_data.csv"
    
    X, y = load_data(file_path)
    
    if X is not None and y is not None:
        standardized_data, scaler = standardize_data(X)
        principal_components, pca, explained_variance_ratio = apply_pca(standardized_data)
        
        print("\nExplained variance ratio of the first two principal components:")
        print(f"PC1: {explained_variance_ratio[0]:.4f} ({explained_variance_ratio[0]*100:.2f}%)")
        print(f"PC2: {explained_variance_ratio[1]:.4f} ({explained_variance_ratio[1]*100:.2f}%)")
        print(f"Total: {sum(explained_variance_ratio):.4f} ({sum(explained_variance_ratio)*100:.2f}%)")
        
        plot_pca_results(principal_components, y, explained_variance_ratio)
        print("\nPCA plot has been saved to 'figures/pca_plot.png'")

if __name__ == "__main__":
    main()