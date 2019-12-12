# Data pre- and post-processing

## Centering and scaling

Use the function `center()` in order to center the data set:

```matlab
[centered_data, ave, centering_name_str] = center(uncentered_data, cent_crit, user_supplied_centering)
```

Use the function `scale()` in order to scale the data set:

```matlab
[scaled_data, scalings, scaling_name_str] = scale(unscaled_data, uncentered_data, scal_crit, user_supplied_scaling)
```

Use the function `uncenter()` in order to uncenter the centered data set:

```matlab
[uncentered_data] = uncenter(centered_data, centerings)
```

Use the function `unscale()` in order to unscale the scaled data set:

```matlab
[unscaled_data] = unscale(scaled_data, scalings)
```
