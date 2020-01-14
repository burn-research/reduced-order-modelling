# ![Logo](documentation/burn_logo.png?thumbnail) Reduced-Order Modelling for combustion data sets

This is a collection of Matlab tools for performing a general Reduced-Order Modelling on data sets. Many of the functions presented here can be used universally on data sets coming from various disciplines but the main focus was to apply the techniques to combustion data sets and hence some methods will be combustion-specific.

The general methodology for the usage of functions in this repository is presented below:

![Screenshot](documentation/rom-methodology.png)

#### General notions that apply across the functions

Whenever we refer to a raw data set we mean a data set that is uncentered and unscaled.

The original raw data set `X_raw` has size `n_obs` x `n_vars`, where `n_obs` is the number of observations and `n_vars` is the number of variables. Typically `n_obs` >> `n_vars` and so `n_vars` determines the dimensionality of the problem.

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/documentation/data-set-for-rom.png">
</p>

This is the shape of the data set that is the input for all the functions presented here. In case your data set is of shape `n_vars` x `n_obs` you need to remember to take its transpose.

# Functions

## Data generation

`Burke_Schumann_solution_CH4_air()`

`add_adaptive_noise_to_data()`

`add_noise_to_data()`

## Data pre-processing

#### Centering and scaling

`center()`

`scale()`

`uncenter()`

`unscale()`

## Clustering

#### Finding clusters

`idx_distance_metrics_clustering()`

`idx_feature_assisted_clustering()`

`idx_kmeans()`

`idx_mixture_fraction_bins()`

`idx_vector_quantization_pca()`

#### Auxiliary functions

`degrade_clusters()`

`get_centroids()`

`get_cluster_populations()`

`get_clusters()`

`get_partition()`

## Reduced-Order Modelling techniques

#### Feature Assisted Clustering (FAC)

#### Independent Component Analysis (ICA)

#### Kernel methods

#### Kriging

#### Local PCA (LPCA)

`lpca()`

`recover_from_lpca()`

#### Non-negative Matrix Factorization (NMF)

#### Polynomial Chaos Expansion (PCE)

## Data post-processing

#### Goodness-of-fit measurements

`average_correlation()`

`cluster_homogeneity_metrics()`

`quality_of_reconstruction_measures()`

#### Data operations

`orthogonalize_factors()`

`trim_factor_weight()`
