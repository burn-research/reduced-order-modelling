# Reduced-Order Modelling techniques

The main methods that are implemented in this repository are presented below.

The techniques implemented here can be applied *globally* (on the entire data set at once) or *locally* (on portions of the entire data set). In order to use the local variants, data set has to be first partitioned into local clusters, for instance using any of the [clustering](https://github.com/burn-research/reduced-order-modelling/tree/master/clustering) techniques available in this repository.

## Principal Component Analysis

[**Global PCA**](https://github.com/burn-research/reduced-order-modelling/tree/master/rom-techniques/pca)

The standard PCA on a pre-processed data set `X` can be run using the built-in Matlab function [`pca()`](https://nl.mathworks.com/help/stats/pca.html):

```matlab
[eigenvectors, scores, eigenvalues, tsquared, variance_explained, mu] = pca(X, 'Centered', false);
```

[**Local PCA**](https://github.com/burn-research/reduced-order-modelling/tree/master/rom-techniques/local-pca)

The local variant of the PCA technique can be found here.

## Non-negative Matrix Factorization

[**Global NMF**](https://github.com/burn-research/reduced-order-modelling/tree/master/rom-techniques/nmf)

The standard PCA on a pre-processed data set `X` can be run using the built-in Matlab function [`nnmf()`](https://nl.mathworks.com/help/stats/nnmf.html):

```matlab
[W,H] = nnmf(X,k);
```

[**Local NMF**](https://github.com/burn-research/reduced-order-modelling/tree/master/rom-techniques/local-nmf)

The local variant of the NMF technique can be found here.

## Proper Orthogonal Decomposition



## Dynamic Mode Decomposition
