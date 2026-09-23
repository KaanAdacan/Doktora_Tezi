# Stage 04 matched 10-ps PCA/clustering policy

PCA and clustering are recomputed from the same 10,000 matched production frames used by Stage 03. MAPK stride = 1 and ADK stride = 5. PCA uses exact covariance eigendecomposition. KMeans k=2..6 uses all 10,000 matched frames. A preferred k is reported only when max Calinski-Harabasz and min Davies-Bouldin agree; otherwise no k or representative is forced.
