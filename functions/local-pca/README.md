# Local Principal Component Analysis (LPCA)

## Local PCA step-by-step

Once you have a division of the original dataset `X` into local clusters saved in a form of a vector `idx` (which can be a result of various clustering techniques such as Vector Quantization (VQ) or K-Means), you can follow these steps to perform PCA in the local clusters.

### Divide the original data set observations to clusters according to `idx`

```matlab
[clusters] = get_clusters(X, idx)
```

### Find PCs, PC-scores and eigenvalues in each cluster

This function performs PCA in each of the local clusters listed in `clusters`.

```matlab
[eigenvectors, scores, eigenvalues, centroids, scales] = lpca(clusters, cent_crit, scal_crit)
```

### Recover original data from low-dimensional representations in each cluster

```matlab
recovered_data = recoverLpca()
```
