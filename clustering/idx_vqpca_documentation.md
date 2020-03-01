# Explanation of the VQPCA algorithm

## Centroids initialization

### Uniform initialization

By default, if no `idx_0` parameter is passed to the function, a uniform centroids initialization will be performed. To explain the elements of the uniform initialization, let's assume that clustering to `k=4` clusters was requested:

<p align="center">
  <img src="https://github.com/burn-research/reduced-order-modelling/raw/master/clustering/dwgs/explanation_of_C_int_and_C.png" width="800">
</p>

The uniform `idx_0` initialization thus assumes that the centroid for every variable is simply that variable's value at the index selected uniformly from the data set.

### Initialization from user-supplied `idx_0`

The user can supply the `idx_0` as one of the parameters and the centroids will be initialized for each cluster as the mean of each variable in that cluster. The function [`get_centroids`](https://github.com/burn-research/reduced-order-modelling/blob/master/clustering/get_centroids.m) is used in that case to obtain the initial cluster centroids.

## Loop

After initializing all necessary parameter the function loops until convergence is reached or until the number of iterations passes the `iter_max` threshold.

Explanation of variables computed inside the loop:

`sq_rec_err` - contains reconstruction error of each observation as if it was assigned to a particular cluster, size `[n_obs, k]`

`C_mat` - for each `j=1:1:k` contains a repeated centroid to form a matrix size `[n_obs, n_vars]`

`[rec_err_min, idx] = min(sq_rec_err, [], 2);` - finds the value and the column-index of the minimum reconstruction error. The index `idx` tells us which cluster the observation should be assigned to to have the minimum reconstruction error.
