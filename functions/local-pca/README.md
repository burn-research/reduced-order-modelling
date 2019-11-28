# Local Principal Component Analysis (LPCA)

## Find division to clusters from LPCA

```matlab
idx = localPCA()
```

## Divide original data set observations to clusters according to `idx`

```matlab
[clusters] = get_clusters(X, idx)
```

## Find PCs, PC-scores and eigenvalues in each cluster

This function performs PCA in each cluster listed in `clusters`.

```matlab
[eigvec, u_scores, eigenvalues, centroids, gamma, eps_rec] = lpca(clusters, cent_crit, scal_crit, is_cpca, idx, cpca_options)
```
