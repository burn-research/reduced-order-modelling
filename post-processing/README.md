# Data post-processing

## Estimation of reconstruction errors

Use the function `quality_of_reconstruction_measures()` to compute reconstruction errors between the original data `X` and the model fit to your data `F`.

```matlab
[r2, nrmse, mae, rmse] = quality_of_reconstruction_measures(X, F)
```

## Data operations

#### Rotation of factors

Factors (eigenvectors or modes) coming from ROM technique can often be rotated to aid the interpretation.

#### Procrustes analysis
