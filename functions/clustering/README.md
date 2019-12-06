# Data clustering

#### Functions starting with the `idx_` prefix

This is a collection of functions that partition (cluster) the data set according to various methods. Each returns a single output: a vector `idx`. Each entry in `idx` specifies the number of cluster to which the observation at that index belongs. An example of such vector could be:

```matlab
idx = [1,1,1,1,2,2,2,2,2,3,3,3,3,3]
```

The above vector has the following interpretation: the first four observations belong to cluster 1, the following five observations belong to cluster 2 and the last five observations belong to cluster 3.

#### Helper functions

This is a collection of functions that do auxiliary tasks in clustering.

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
