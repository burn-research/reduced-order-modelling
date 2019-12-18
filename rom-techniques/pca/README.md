# Principal Component Analysis (PCA)

The standard PCA on a pre-processed data set `X` can be run using the built-in Matlab function [`pca()`](https://nl.mathworks.com/help/stats/pca.html):

```matlab
[eigenvectors, scores, eigenvalues, tsquared, variance_explained, mu] = pca(X, 'Centered', false);
```

## Functions descriptions

### Recover the original data set

Use the function `recover_from_pca()` to recover a rank-`q` approximation to the original raw data set `X_raw`.

```matlab
[X_app_pca] = recover_from_pca(eigenvectors, scores, q, centerings, scalings)
```

If `q` is not specified, the highest rank approximation will be made. If `centerings` or `scalings` is not specified, the approximated data set will not be uncentered or unscaled at the end.

![Screenshot](https://raw.githubusercontent.com/camillejr/ulb-atm-phd/master/PCA/DWGs/PCA-example-subplot.png)
