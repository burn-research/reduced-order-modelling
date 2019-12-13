# Data post-processing

## Estimation of reconstruction errors

Use the function `quality_of_reconstruction_measures()` to compute reconstruction errors between the original data `X` and the model fit to your data `F`.

```matlab
[r2, nrmse, mae, rmse] = quality_of_reconstruction_measures(X, F)
```

## Data operations

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

#### Procrustes analysis

The Procrustes analysis can be used to compare two low-dimensional manifolds (e.g. coming from two ROM techniques). The analysis can be performed using the built-in Matlab function [`procrustes()`](https://nl.mathworks.com/help/stats/procrustes.html):

```Matlab
[d, Z, transform] = procrustes(X,Y)
```

where `X` an `Y` are manifolds that need to be conformed.

#### Cluster homogeneity metrics

In order to measure the homogeneity inside clusters found from various clustering techniques, use the function `cluster_homogeneity_metrics()`:

```Matlab
[mean_rmse, mean_silhouette, mean_delta, var_delta] = cluster_homogeneity_metrics(X, idx, n_pcs, cent_crit, scal_crit)
```
