# Local Principal Component Analysis (LPCA)


## Local PCA step-by-step

Once you have a division of a dataset into local clusters saved in a form of a vector `idx`, you can follow these steps to perform PCA in the local clusters.

### Divide the original data set observations to clusters according to `idx`

```matlab
[clusters] = get_clusters(X, idx)
```

### Find PCs, PC-scores and eigenvalues in each cluster

This function performs PCA in each cluster listed in `clusters`.

```matlab
[eigvec, u_scores, eigenvalues, centroids, scales] = lpca(clusters, cent_crit, scal_crit)
```

### Recover original data from low-dimensional representations in each cluster

```matlab
recovered_data = recoverLpca()
```
