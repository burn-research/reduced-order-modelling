# ![Logo](documentation/burn_logo.png?thumbnail) Reduced-Order Modelling for combustion data sets

This is a collection of Matlab tools for performing a general Reduced-Order Modelling on data sets.

Many of the functions presented here can be used universally on data sets coming from various disciplines but the main focus was to apply the techniques to combustion data sets and hence some methods are combustion-specific.

The general workflow for the usage of functions in this repository is presented below:

![Screenshot](documentation/rom-methodology.png)

> **Figure 1.**

# Functions

Each of the boxes present in **Figure 1.** contains its functions for performing specific tasks. You can explore the functions and documentation for each of the boxes by clicking the relevant section below:

### [`Data generation`](https://github.com/burn-research/reduced-order-modelling/tree/master/data-generation)

### [`Data pre-processing`](https://github.com/burn-research/reduced-order-modelling/tree/master/pre-processing)

### [`Clustering`](https://github.com/burn-research/reduced-order-modelling/tree/master/clustering)

### [`ROM techniques`](https://github.com/burn-research/reduced-order-modelling/tree/master/rom-techniques)

### [`Data post-processing`](https://github.com/burn-research/reduced-order-modelling/tree/master/post-processing)

## General notions that apply across the functions in this repository

Whenever we refer to a raw data set we mean a data set that is uncentered and unscaled.

The raw data set `X_raw` has size `n_obs` x `n_vars`, where `n_obs` is the number of observations and `n_vars` is the number of variables. Typically `n_obs` >> `n_vars` and so `n_vars` determines the dimensionality of the problem.

<p align="center">
<img src="https://github.com/burn-research/reduced-order-modelling/raw/master/documentation/data-set-for-rom.png">
</p>

> **Figure 2.**

This is the shape of the data set that is the input for all the functions presented here. In case your data set is of shape `n_vars` x `n_obs` you need to remember to take its transpose.
