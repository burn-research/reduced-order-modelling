# Local Principal Component Analysis (LPCA)

The figure below presents the difference between applying PCA globally vs. locally.

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/documentation/global-vs-local-pca-subplot.png" width="300">
</p>

## Local PCA step-by-step

This is the typical workflow for performing Local PCA.

### 1. Obtain the data set partitioning `idx`

To cluster the data set you can use any technique that you want. The ultimate goal of clustering is to obtain the vector `idx` which classifies every observations to a particular cluster. In this repository few techniques are implemented in the [`clustering`](https://github.com/burn-research/reduced-order-modelling/tree/master/clustering) directory:

- Vector Quantization PCA (VQPCA)
- Mixture fraction clustering
- K-Means clustering
- Feature Assisted Clustering (FAC)

### 2. Pre-process the data set (center and scale)

The LPCA typically starts with pre-processing the raw data set `X_raw` (uncentered and unscaled).

```matlab
[X, centerings] = center(X_raw, cent_crit);
[X, scalings] = scale(X, X_raw, scal_crit);
```

The same `centerings` and `scalings` will be used to uncenter and unscale the reconstructed data set within the function `recover_from_lpca()`.

### 3. Divide the original data set observations to clusters according to `idx`

Once you have a classification of the data set `X` into clusters you can divide the data set observations into clusters according to `idx`:

```matlab
[clusters] = get_clusters(X, idx)
```

### 4. Find PCs, PC-scores and eigenvalues in each cluster

This function performs Principal Component Analysis in each of the local clusters listed in `clusters`:

```matlab
[eigenvectors, scores, eigenvalues, centroids, local_scales] = lpca(clusters, q, cent_crit, scal_crit)
```

The parameter `q` specifies how many components should be kept within the returned variables. If `q` is not specified, all components are returned (and the total number of components is the number of variables in a data set).

Data samples within each cluster are centered with a selected `cent_crit` criterion. If `cent_crit` is not specified, mean-centering is performed by default. The centerings are returned in the `centroids` variable.

Data samples within each cluster can also be scaled with a selected `scal_crit` criterion. If `scal_crit` is not specified, data samples are not scaled within each cluster. The scalings are returned in the `local_scales` variable.

If defaults are used for all possible variables, the function can only be run with the `clusters` variable:

```matlab
[eigenvectors, scores, eigenvalues, centroids, local_scales] = lpca(clusters)
```

### 5. Recover the original data from low-dimensional representations in each cluster

In order to reconstruct the original data set `X_raw`, the following function can be used:

```matlab
[X_app_lpca] = recover_from_lpca(idx, eigenvectors, scores, q, centroids, local_scalings, centerings, scalings)
```
