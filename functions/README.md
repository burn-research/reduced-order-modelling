# Reduced-Order Modelling functions

This is a library of Matlab functions for performing Reduced-Order Modelling (ROM).

The general methodology for the usage of functions in this repository is presented on the graph below:

![Screenshot](rom-methodology.png)

## Available techniques

### Data pre- and post-processing

#### Centering and scaling

`center()`

`scale()`

`uncenter()`

`unscale()`

### Clustering

`get_centroids()`

`get_clusters()`

`get_cluster_populations()`

`idx_mixture_fraction_bins()`

### Principal Component Analysis (PCA)

#### Local PCA (LPCA, VQPCA)

`localPCA()`

`lpca()`

#### Variable PCA (FPCA)

`FPCA()`

#### Kernel PCA (KPCA)

### Independent Component Analysis (ICA)

`fastICA()`

`localICA()`

### Non-negative Matrix Factorization (NMF)

`localNNMF()`

### Kriging

`performKriging()`

### Quality of reconstruction measurements

`r_square()`
