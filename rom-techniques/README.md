# Reduced-Order Modelling techniques

The main methods that are implemented in this repository are presented below.

## Principal Component Analysis

**Global PCA**

The standard PCA on a pre-processed data set `X` can be run using the built-in Matlab function [`pca()`](https://nl.mathworks.com/help/stats/pca.html):

```matlab
[eigenvectors, scores, eigenvalues, tsquared, variance_explained, mu] = pca(X, 'Centered', false);
```

**Local PCA**

PCA can also be applied in local clusters that can be found by various [clustering](https://github.com/burn-research/reduced-order-modelling/tree/master/clustering) techniques. Some of the commonly used partitioning algorithms in combustion data sets are: FPCA, where data is partitioned based on mixture fraction bins or VQPCA, where data is partitioned according to the Vector Quantization (VQ) algorithm. Once the partitioning is found, the local PCA can be applied and the functions are available in the [local-pca](https://github.com/burn-research/reduced-order-modelling/tree/master/rom-techniques/local-pca) directory.

## Independent Component Analysis


## Non-negative Matrix Factorization


## Kernel methods


## Kriging


## Polynomial Chaos Expansion
