# Utility functions

## `CoupledSDEs` utility functions

```@docs
covariance_matrix
diffusion_matrix
dynamic_rule
noise_process
current_state
set_state!
drift
div_drift
```

### `length(sys.u)`-dimensional box

```@docs
intervals_to_box(bmin::Vector, bmax::Vector)
```