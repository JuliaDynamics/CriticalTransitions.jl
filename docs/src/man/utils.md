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

## Saving data

```@docs
make_jld2(text::String, relpath::String="")
make_h5(text::String, relpath::String="")
```

### `length(sys.u)`-dimensional box

```@docs
intervals_to_box(bmin::Vector, bmax::Vector)
```