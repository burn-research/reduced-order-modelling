# Data clustering

#### Functions starting with the `idx_` prefix

This is a collection of functions that partition (cluster) the data set according to various methods. Each returns a single output: a vector `idx`. Each entry in `idx` specifies the number of cluster to which the observation at that index belongs. An example of such vector could be:

```matlab
idx = [1,1,1,1,2,2,2,2,2,3,3,3,3,3]
```

The above vector has the following interpretation: the first four observations belong to cluster 1, the following five observations belong to cluster 2 and the last five observations belong to cluster 3.

## Function descriptions

### Bins of mixture fraction

Use the function `idx_mixture_fraction_bins()` to partition the data into `k` clusters (bins) according to the mixture fraction space `Z`. Bins intervals are found starting from the first one which divides the `Z` space into two parts at the stoichiometric mixture fraction `Z_stoich`.

```matlab
[idx] = idx_mixture_fraction_bins(Z, k, Z_stoich)
```

#### Related functions

In order to plot the visualisation of this clustering technique, use the function [`plot_mixture_fraction_divided_to_clusters`](https://github.com/burn-research/plotting/blob/master/plot_mixture_fraction/plot_mixture_fraction_divided_to_clusters.m) from the [`plotting`](https://github.com/burn-research/plotting) repository:

```matlab
plot_mixture_fraction_divided_to_clusters(data, Z, idx, np, destination)
```

A visualisation of the partitioning according to this method:

![Screenshot](dwgs/idx_mixture_fraction_bins.png)

### K-Means

# Helper functions

This is a collection of functions that do auxiliary tasks in clustering.

## Function descriptions

### `get_centroids`

This functions returns a matrix of centroids of every cluster. Centroids are computed as the mean of all the observations of a specific variable in a particular cluster.

```matlab
[centroids] = get_centroids(X, idx)
```

### `get_clusters`

This functions returns a cell array of data set divided into clusters.

```matlab
[clusters] = get_clusters(X, idx)
```

### `get_cluster_populations`

This functions returns a vector whose each entry is a number of observations in each cluster.

```matlab
[populations] = get_cluster_populations(idx)
```

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
