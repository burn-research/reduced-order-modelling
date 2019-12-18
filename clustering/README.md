# Data clustering

#### Functions starting with the `idx_` prefix

This is a collection of functions that partition (cluster) the data set according to various techniques. Each function returns a single output: an integer vector `idx`. Each entry inside `idx` specifies the number of cluster to which the observation at that index was classified to. An example of such vector could be:

```matlab
idx = [1,1,1,1,2,2,2,2,2,3,3]
```

The above vector has the following interpretation: the first four observations belong to cluster 1, the following five observations belong to cluster 2 and the last two observations belong to cluster 3.

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/documentation/idx-X.png" width="300">
</p>

## Clustering functions descriptions

### Bins of mixture fraction

Use the function `idx_mixture_fraction_bins()` to partition the data into `k` clusters (bins) according to the mixture fraction space `Z`. Bins intervals are found starting from the first split which divides the `Z` space into two parts at the stoichiometric mixture fraction `Z_stoich`. Next, the **rich** and **lean** sides are divided equally into clusters.

```matlab
[idx] = idx_mixture_fraction_bins(Z, k, Z_stoich)
```

An example visualisation of the partitioning according to this method:

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/clustering/dwgs/idx_mixture_fraction_bins.png" width="500">
</p>

### Vector Quantization Principal Component Analysis (VQPCA)

Use the function `idx_vector_quantization_pca()` to partition the data into `k` clusters according to the Vector Quantization (VQ) algorithm.

```matlab
[idx] = idx_vector_quantization_pca(X, n_eigs, k, cent_crit, scal_crit, idx_0)
```

### Feature Assisted Clustering (FAC)

Use the function `idx_feature_assisted_clustering()` to partition the data into `k` clusters according to the Feature Assisted Clustering.

```matlab
idx] = idx_feature_assisted_clustering(X, artificial_modes, cent_crit, scal_crit)
```

### K-Means

Use the function `idx_kmeans()` to partition the data into `k` clusters according to the K-Means algorithm.

```matlab
[idx] = idx_kmeans(X, k)
```

> ### Related functions
>
> In order to plot a visualisation of a clustering technique, use the function [`plot_mixture_fraction_divided_to_clusters`](https://github.com/burn-research/plotting/blob/master/plot_mixture_fraction/plot_mixture_fraction_divided_to_clusters.m) from the [`plotting`](https://github.com/burn-research/plotting) repository:
>
> ```matlab
plot_mixture_fraction_divided_to_clusters(data, Z, idx, np, destination)
```

# Auxiliary functions

This is a collection of functions that do auxiliary tasks in clustering.

## Function descriptions

### `degrade_clusters`

This function degrades cluster numeration to consecutive integers. For example, if the initial `idx` is:

```matlab
idx = [1,1,1,2,2,4,4,4,2]
```

the new `idx` returned by the `degrade_clusters()` function will be:

```matlab
idx = [1,1,1,2,2,3,3,3,2]
```

so that the cluster numeration is composed of consecutive integers.

```matlab
[new_idx] = degrade_clusters(idx)
```

### `get_centroids`

This functions returns a matrix of centroids of every cluster. Centroids are computed as the mean of all the observations of a specific variable in a particular cluster.

```matlab
[centroids] = get_centroids(X, idx)
```

### `get_cluster_populations`

This functions returns a vector whose each entry is a number of observations in each cluster.

```matlab
[populations] = get_cluster_populations(idx)
```

### `get_clusters`

This functions returns a cell array of data set divided into clusters.

```matlab
[clusters] = get_clusters(X, idx)
```

### `get_partition`

This function partitions the dataset to clusters given by `idx` vector. It returns a cell array of data set divided into clusters and a cell array of indices of original observations diveded into clusters. If less than `n_vars` observations are assigned to a particular cluster, this cluster will be removed.

```matlab
[clusters, clusters_idx, k_new] = get_partition(X, idx, k)
```
