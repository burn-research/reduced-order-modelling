# Local Principal Component Analysis (LPCA)

## Local PCA step-by-step

This is the typical workflow when performing LPCA.

### Obtain the data set partitioning `idx`

To cluster the data set you can use any technique that you want. The ultimate goal of clustering is to obtain the vector `idx` which classifies every observations to a particular cluster. In this repository three techniques are implemented in the [`clustering`](https://github.com/burn-research/reduced-order-modelling/tree/master/clustering) directory:

- Vector Quantization PCA (VQPCA)
- K-Means clustering
- Mixture fraction clustering

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/documentation/idx-X.png" width="300">
</p>

### Pre-process the data set (center and scale)

The LPCA typically starts with pre-processing the global raw data set `X_raw`.

```matlab
[X, centerings] = center(X_raw, cent_crit);
[X, scalings] = scale(X, X_raw, scal_crit);
```

The same `centerings` and `scalings` will later be used to uncenter and unscale the reconstructed data set within the function `recover_from_lpca()`.

### Divide the original data set observations to clusters according to `idx`

Once you have a division of the original dataset `X` into local clusters saved in a form of a vector `idx` you can divide the data set observations into clusters according to `idx`:

```matlab
[clusters] = get_clusters(X, idx)
```

### Find PCs, PC-scores and eigenvalues in each cluster

This function performs PCA in each of the local clusters listed in `clusters`.

```matlab
[eigenvectors, scores, eigenvalues, centroids, local_scales] = lpca(clusters, q, cent_crit, scal_crit)
```

The parameter `q` specifies how many components should be kept within the returned variables. If `q` is not specified, all components are returned (the total number of components is the number of variables in a data set).

Data samples within each cluster are centered with a selected `cent_crit` criterion. If `cent_crit` is not specified, mean-centering is performed by default. The centerings are returned in the `centroids` variable.

Data samples within each cluster can also be scaled with a selected `scal_crit` criterion. If `scal_crit` is not specified, data samples are not scaled within each cluster. The scalings are returned in the `local_scales` variable.

If defaults are used for all variables, the function can only be run with the `clusters` variable:

```matlab
[eigenvectors, scores, eigenvalues, centroids, local_scales] = lpca(clusters)
```

### Recover the original data from low-dimensional representations in each cluster

In order to reconstruct

```matlab
[X_app_lpca] = recover_from_lpca(idx, eigenvectors, scores, q, centroids, local_scalings, centerings, scalings)
```
