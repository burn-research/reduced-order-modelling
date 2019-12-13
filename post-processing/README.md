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

#### Procrustes analysis

The Procrustes analysis can be used to compare two low-dimensional manifolds (e.g. coming from two ROM techniques). The analysis can be performed using the built-in Matlab function [`procrustes()`](https://nl.mathworks.com/help/stats/procrustes.html):

```Matlab
[d, Z, transform] = procrustes(X,Y)
```

where `X` an `Y` are manifolds that need to be conformed.
