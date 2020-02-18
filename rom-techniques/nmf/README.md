# Non-negative Matrix Factorization (NMF)

The standard NMF on a pre-processed data set `X` can be run using the built-in Matlab function [`nnmf()`](https://nl.mathworks.com/help/stats/nnmf.html):

```matlab
[W, H] = nnmf(X, q);
```
