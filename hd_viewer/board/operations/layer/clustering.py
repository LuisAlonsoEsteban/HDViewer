
import numpy as np

from sklearn.cluster import KMeans, DBSCAN, HDBSCAN, AgglomerativeClustering
from sklearn.mixture import GaussianMixture

ALGORITHMS = ("leiden", "kmeans", "dbscan", "hdbscan", "gaussian_mixture", "agglomerative")

def fit_clusters(values, algorithm, *, n_clusters=8, eps=0.5, min_samples=5,
                 min_cluster_size=5, covariance_type="full", linkage="ward", random_state=0):
                 
    if algorithm not in ALGORITHMS or algorithm == "leiden":
        raise ValueError(f"Unsupported PCA clustering algorithm: {algorithm}")

    if not np.isfinite(values).all():
        raise ValueError("Clustering requires finite PCA coordinates.")

    if algorithm in ("kmeans", "gaussian_mixture", "agglomerative"):
        if not isinstance(n_clusters, (int, np.integer)) or not 1 <= n_clusters <= len(values):
            raise ValueError("Number of clusters must be between 1 and the number of observations.")
    
    estimators = {
        "kmeans": lambda: KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10),
        "dbscan": lambda: DBSCAN(eps=eps, min_samples=min_samples),
        "hdbscan": lambda: HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples),
        "gaussian_mixture": lambda: GaussianMixture(n_components=n_clusters,
            covariance_type=covariance_type, random_state=random_state),
        "agglomerative": lambda: AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage),
    }
    
    return estimators[algorithm]().fit_predict(values)
