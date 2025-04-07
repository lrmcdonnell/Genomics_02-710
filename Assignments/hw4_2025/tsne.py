import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)

def load_data(file_path):
    try:
        data = pd.read_csv(file_path)
        print(f"Data loaded successfully with shape: {data.shape}")
        return data
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def standardize_data(data):
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(data)
    return standardized_data, scaler

def apply_tsne(standardized_data, perplexity_values, learning_rate=200):
    tsne_results = {}
    
    for perplexity in perplexity_values:
        print(f"Applying t-SNE with perplexity={perplexity}, learning_rate={learning_rate}")
        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=1000,
            random_state=42
        )
        
        tsne_result = tsne.fit_transform(standardized_data)
        tsne_results[perplexity] = tsne_result
    
    return tsne_results

def plot_tsne_results(tsne_results, labels, perplexity_values):
    fig, axes = plt.subplots(1, len(perplexity_values), figsize=(18, 6))
    
    for i, perplexity in enumerate(perplexity_values):
        ax = axes[i]
        tsne_result = tsne_results[perplexity]
        
        scatter = ax.scatter(
            tsne_result[:, 0],
            tsne_result[:, 1],
            c=labels,
            cmap='viridis',
            alpha=0.8,
            s=50
        )
        
        legend1 = ax.legend(*scatter.legend_elements(),
                            title="Classes", loc="upper right")
        ax.add_artist(legend1)
        
        ax.set_title(f't-SNE with Perplexity={perplexity}')
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('tsne_results.png', dpi=300, bbox_inches='tight')
    plt.show()

def calculate_silhouette_scores(tsne_results, labels, perplexity_values):
    silhouette_scores = {}
    
    for perplexity in perplexity_values:
        tsne_result = tsne_results[perplexity]
        score = silhouette_score(tsne_result, labels)
        silhouette_scores[perplexity] = score
        print(f"Silhouette score for perplexity={perplexity}: {score:.4f}")
    
    best_perplexity = max(silhouette_scores, key=silhouette_scores.get)
    print("\n=== Silhouette Score Results ===")
    print("Perplexity\tSilhouette Score")
    print("-----------------------------------")
    for perplexity in sorted(perplexity_values):
        score = silhouette_scores[perplexity]
        is_best = "  ★ BEST" if perplexity == best_perplexity else ""
        print(f"{perplexity}\t\t{score:.4f}{is_best}")
    print("-----------------------------------")
    print(f"Best perplexity value: {best_perplexity} (Silhouette score: {silhouette_scores[best_perplexity]:.4f})")
    
    return silhouette_scores, best_perplexity

def main():
    file_path = "provided_data/gene_expression_data.csv"
    
    data = load_data(file_path)
    
    if data is not None:
        print("First 5 rows of data:")
        print(data.head())
        
        potential_label_cols = []
        for col in data.columns:
            if col.lower() in ['label', 'class', 'type', 'category']:
                potential_label_cols.append(col)
        
        gene_data = data
        
        if potential_label_cols:
            label_col = potential_label_cols[0]
            print(f"Found label column: {label_col}")
            labels = data[label_col].values
            gene_data = data.drop(label_col, axis=1)
        
        elif data.iloc[:, 0].nunique() <= 5:
            print("First column appears to contain labels")
            labels = data.iloc[:, 0].values
            gene_data = data.iloc[:, 1:]
        
        elif any(isinstance(x, str) for x in data.iloc[0].values):
            print("First row appears to contain labels")
            
            labels = data.iloc[0, 1:].values
            gene_data = data.iloc[1:, 1:].reset_index(drop=True)
            gene_data = gene_data.astype(float)
            
            unique_labels = np.unique([str(x) for x in labels])
            label_map = {label: i for i, label in enumerate(unique_labels)}
            labels_int = np.array([label_map.get(str(x), 0) for x in labels])
        
        else:
            print("Could not detect labels, generating dummy labels")
            n_samples = data.shape[0]
            labels_int = np.concatenate([
                np.zeros(n_samples // 2, dtype=int),
                np.ones(n_samples - n_samples // 2, dtype=int)
            ])
            gene_data = data
        
        if 'labels_int' not in locals():
            if isinstance(labels[0], (str, bool)):
                unique_labels = np.unique(labels)
                label_map = {label: i for i, label in enumerate(unique_labels)}
                labels_int = np.array([label_map.get(label, 0) for label in labels])
            else:
                labels_int = labels.astype(int)
        
        print(f"Number of samples: {len(labels_int)}")
        print(f"Label distribution: {np.bincount(labels_int)}")
        
        print(f"Gene data shape before standardization: {gene_data.shape}")
        standardized_data, _ = standardize_data(gene_data)
        print(f"Standardized data shape: {standardized_data.shape}")
        
        perplexity_values = [5, 10, 50]
        
        tsne_results = apply_tsne(standardized_data, perplexity_values)
        
        plot_tsne_results(tsne_results, labels_int, perplexity_values)
        
        silhouette_scores, best_perplexity = calculate_silhouette_scores(
            tsne_results, labels_int, perplexity_values
        )
        
        plt.figure(figsize=(10, 6))
        plt.bar(
            [str(p) for p in perplexity_values],
            [silhouette_scores[p] for p in perplexity_values],
            color='skyblue'
        )
        plt.xlabel('Perplexity')
        plt.ylabel('Silhouette Score')
        plt.title('Silhouette Scores for Different Perplexity Values')
        plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig('silhouette_scores.png', dpi=300, bbox_inches='tight')
        plt.show()

if __name__ == "__main__":
    main()