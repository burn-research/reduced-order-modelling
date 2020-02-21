# Explanation of the elements in the VQPCA algorithm

## Centroids initialisation

### Uniform initialisation

To explain the elements of the uniform initialisation, let's assume that clustering to `k=4` clusters was requested.

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/clustering/dwgs/explanation_of_C_int_and_C.png" width="800">
</p>

## Loop

`sq_rec_err`

`C_mat` - for each `j=1:1:k` contains a repeated centroid to form a matrix size `[n_obs, n_vars]`

`sq_rec_err` - contains reconstruction error of each observation if it was assigned to a particular cluster, size `[n_obs, k]`

`[rec_err_min, idx] = min(sq_rec_err, [], 2);` - finds the value and the column-index of the minimum reconstruction error. The index `idx` tells us which cluster the observation should be assigned to to have the minimum reconstruction error.
