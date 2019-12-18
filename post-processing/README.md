# Data post-processing

## Metrics

### Reconstruction quality

Use the function `quality_of_reconstruction_measures()` to compute reconstruction errors between the original data `X` and the model fit to your data `F`.

```matlab
[r2, nrmse, mae, rmse] = quality_of_reconstruction_measures(X, F)
```

### Clustering quality

#### Cluster homogeneity metrics

In order to measure the homogeneity inside clusters found from various clustering techniques, use the function `cluster_homogeneity_metrics()`:

```Matlab
[mean_rmse, mean_silhouette, mean_delta, var_delta] = cluster_homogeneity_metrics(X, idx, q, cent_crit, scal_crit)
```

## Data operations

### Operations on factors/modes/eigenvectors

#### Rotation of factors

Factors (eigenvectors or modes) coming from ROM technique can be rotated to aid the interpretation.

This can be done using the built-in Matlab function [`rotatefactors()`](https://nl.mathworks.com/help/stats/rotatefactors.html):

```Matlab
B = rotatefactors(A)
```

One of the possible methods is a Varimax rotation:

```Matlab
B = rotatefactors(A, 'Method', 'varimax')
```

#### Factor weight trimming

In order to further aid the interpretation of factors, small weights on factors can be set to zero. Use the function `trim_factor_weight()` to set the weights below a `threshold` to zero, or the function `trim_factor_weight_iteratively()`, which trims weights iteratively and re-normalizes the factor after each trim. This function thus will never return a zero vector, which might happen with `trim_factor_weight` if you use too high threshold.

Use `trim_factor_weight_iteratively()` when your threshold is of the order of the highest element of a vector.

#### Orthogonalize factors

To orthogonalize factors, use the function `orthogonalize_factors()` which performs the Gram-Schmidt orthogonalization.

#### Unique modes

This is a set of functions that help performing analysis of unique/repeating modes between modes find by global techniques and those find in local clusters.

Below is a description of possible outcomes.

In all codes below:

`k` is the number of clusters.

`alpha` is the number between 0 and 1 that specifies the level of similarity between the modes that you desire. The similarity is measured on the basis of a dot (inner) product between modes and alpha is the condition for the cutoff value of the dot products. Example: set to 0.8 if you want 80% similarity between the modes as the cutoff value.

##### Which modes from local clusters are unique and do not repeat within the global modes?

Use the function `unique_modes_clusters()` with the input parameters: `modes_global`, `modes_clusters`, `k` and `alpha`.

The outputs of the code are:

`uni_modes_clusters` is the matrix containing unique modes from local clusters that are not present within the global modes.

`mode_annotation` is a vector whose each entry is associated to a respective column of `uni_modes_clusters` and specifies the number of the mode.

`cluster_annotation` is a vector whose each entry is associated to a respective column of `uni_modes_clusters` and specifies in which cluster the unique mode was present.

##### Example

```matlab
clc, clear, close all

k = 3;
alpha = 0.99;

g1 = [1;2;1;2]/norm([1;2;1;2]);
g2 = [2;3;3;2]/norm([2;3;3;2]);
g3 = [3;5;3;1]/norm([3;5;3;1]);
g4 = [0;4;2;4]/norm([0;4;2;4]);

l1 = [5;5;5;5]/norm([5;5;5;5]);
l2 = [3;4;6;6]/norm([3;4;6;6]);
l3 = [7;0;7;7]/norm([7;0;7;7]);

modes_global = [g1, g2, g3, g4, g2, g3, g4];
modes_local = [l1, l2, l3, g2, g2, l1];

[uni_modes_clusters, cluster_annotation] = unique_modes_clusters(modes_global, modes_clusters, k, alpha);
```

The above code will produce the following results:

`uni_modes_clusters = [l1, l2, l3, l1]` - these are all the unique vectors that make up `modes_clusters`.

`mode_annotation = [1, 2, 1, 2]` - there are 4 vectors that do not repeat within the global modes and they are modes number: 1, 2, 1 and 2 within their own clusters respectively. (`l1` is the **first** mode in the first cluster, `l2` is the **second** mode in the first cluster and so on...)

`clusters_annotation = [1, 1, 2, 3]` - there are 4 vectors that do not repeat within the global modes and they belong to clusters: 1, 1, 2 and 3 respectively. (`l3` is the first mode in the **second** cluster, `l1` is also the second mode in the **third** cluster.)

The code will also print verbally:

> *From 6 local modes, 4 are not repeating within the global modes.*

##### How many and which global modes are repeating in local clusters?

Use the function `repeating_modes_global()` with the input parameters: `modes_global`, `modes_local`, `k` and `alpha`.

The outputs of the code are:

`rep_modes_global` is the matrix containing unique global modes which repeat in local clusters.

`rep_modes_local` is the matrix with all the modes from local clusters which are the same as any of the global modes.

`cluster_annotation` is a vector whose each entry is associated to a respective column of `rep_modes_local` and specifies in which cluster this repeating mode was present.

`no_repeating` is a scalar specifying how many unique repetition from global modes are also present in local clusters. It will always be equal to the number of columns of `rep_modes_global`.

##### Example

```matlab
k = 3;
alpha = 0.95;

g1 = [1;2;1;2]/norm([1;2;1;2]);
g2 = [2;3;3;2]/norm([2;3;3;2]);
g3 = [3;5;3;1]/norm([3;5;3;1]);
g4 = [0;4;2;4]/norm([0;4;2;4]);

l1 = [5;5;5;5]/norm([5;5;5;5]);
l2 = [3;4;6;6]/norm([3;4;6;6]);
l3 = [7;0;7;7]/norm([7;0;7;7]);

modes_global = [g1, g2, g3, g4, g2, g3, g4];
modes_local = [l1, l2, l3, g2, g2, l1];

[rep_modes_global, rep_modes_local, cluster_annotation, no_repeating] = repeating_modes_global(modes_global, modes_local, k, alpha);
```

The above code will produce the following results:

`rep_modes_global = [g2]` - there is 1 unique repeating vector (`g2`) inside global modes which repeats in local modes (**it will be counted only once**, even though it was present twice in `modes_global`).

`rep_modes_local = [g2, g2]` - this 1 unique repeating vector (`g2`) from global modes is present twice in local modes.

`cluster_annotation = [2, 3]` - this 1 unique repeating vector (`g2`) from global modes is present in cluster 2 and cluster 3.

`no_repeating = 1` - there is 1 unique repeating vector (`g2`) inside global modes which repeats in local modes.

##### How many unique modes there are within local clusters?

Use the function `unique_modes_within_clusters()` with the input parameters: `modes_clusters`, `k` and `alpha`.

The outputs of the code are:

`unique_modes` is the matrix containing unique modes that make up all the modes in local clusters.

`no_unique` is a scalar specifying how many unique modes make up all the modes in local clusters.

##### Example

```matlab
k = 2;
alpha = 0.99;

c1 = [1;2;1;2]/norm([1;2;1;2]);
c2 = [2;3;3;2]/norm([2;3;3;2]);
c3 = [3;5;3;1]/norm([3;5;3;1]);
c4 = [0;4;2;4]/norm([0;4;2;4]);
c5 = [5;5;5;5]/norm([5;5;5;5]);
c6 = [3;4;6;6]/norm([3;4;6;6]);
c7 = [7;0;7;7]/norm([7;0;7;7]);

modes_clusters = [c1, c6, c2, c3, c4, c2, c3, c4, c5, c5, c7, c1];

[unique_modes, no_unique] = unique_modes_within_clusters(modes_clusters, k, alpha);
```

The above code will produce the following results:

`unique_modes = [c1, c6, c2, c3, c4, c5, c7]` - these are all the unique vectors that make up `modes_clusters`.

`no_unique = 7` - there are 7 unique vectors (this is always equal to the number of columns of `unique_modes`.

The code will also print verbally:

> *Out of 12 modes in local clusters, there are 7 unique modes.*

### Operations on low-dimensional manifolds

#### Procrustes analysis

The Procrustes analysis can be used to compare two low-dimensional manifolds (e.g. coming from two ROM techniques). The analysis can be performed using the built-in Matlab function [`procrustes()`](https://nl.mathworks.com/help/stats/procrustes.html):

```Matlab
[d, Z, transform] = procrustes(X,Y)
```

where `X` an `Y` are manifolds that need to be conformed.
