struct TransitionPathEnsemble
    paths::Vector
    times::Vector
    success_rate::Real
    t_res::Real
    t_trans::Real
end;

function prettyprint(tpe::TransitionPathEnsemble)
    "Transition path ensemble of $(length(tpe.times)) samples
    - sampling success rate:      $(round(tpe.success_rate, digits=3))
    - mean residence time:        $(round(tpe.t_res, digits=3))
    - mean transition time:       $(round(tpe.t_trans, digits=3))
    - normalized transition rate: $(round(tpe.t_res/tpe.t_trans, digits=1))"
end

Base.show(io::IO, tpe::TransitionPathEnsemble) = print(io, prettyprint(tpe))

# """
#     transition(sys::StochSystem, x_i::State, x_f::State; kwargs...)
# Generates a sample transition from point `x_i` to point `x_f`.

# This function simulates `sys` in time, starting from initial condition `x_i`, until entering a `length(sys.u)`-dimensional ball of radius `rad_f` around `x_f`.

# ## Keyword arguments
# * `rad_i=0.1`: radius of ball around `x_i`
# * `rad_f=0.1`: radius of ball around `x_f`
# * `cut_start=true`: if `false`, returns the whole trajectory up to the transition
# * `dt=0.01`: time step of integration
# * `tmax=1e3`: maximum time when the simulation stops even `x_f` has not been reached
# * `rad_dims=1:length(sys.u)`: the directions in phase space to consider when calculating the radii
#   `rad_i` and `rad_f`. Defaults to all directions. To consider only a subspace of state space,
#   insert a vector of indices of the dimensions to be included.
# * `solver=EM()`: numerical solver. Defaults to Euler-Mayurama.
# * `progress`: shows a progress bar with respect to `tmax`

# ## Output
# `[path, times, success]`
# * `path` (Matrix): transition path (size [dim × N], where N is the number of time points)
# * `times` (Vector): time values (since start of simulation) of the path points (size N)
# * `success` (bool): if `true`, a transition occured (i.e. the ball around `x_f` has been reached), else `false`
# * `kwargs...`: keyword arguments passed to [`simulate`](@ref)

# See also [`transitions`](@ref), [`simulate`](@ref).
# """
function transition(sys::CoupledSDEs, x_i, x_f;
        rad_i = 0.1,
        rad_f = 0.1,
        tmax = 1e3,
        cut_start = true,
        rad_dims = 1:length(current_state(sys)),
        kwargs...)
    condition(u, t, integrator) = subnorm(u - x_f; directions = rad_dims) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)

    sim = simulate(sys, tmax, x_i; callback = cb_ball, kwargs...)
    success = sim.retcode == SciMLBase.ReturnCode.Terminated

    simt = sim.t
    if cut_start
        idx = size(sim)[2]
        dist = norm(sim[:, idx] - x_i)
        while dist > rad_i
            idx -= 1
            dist = norm(sim[:, idx] - x_i)
            idx < 1 &&
                error("Trajactory never left the initial state sphere. Increase tmax or decrease rad_i.")
        end
        sim = sim[:, idx:end]
        simt = simt[idx:end]
    end

    sim, simt, success
end;

"""
    transitions(sys::StochSystem, x_i::State, x_f::State, N=1; kwargs...)

Generates an ensemble of `N` transition samples of `sys` from point `x_i` to point `x_f`.

This function repeatedly calls the [`transition`](@ref) function to efficiently generate an ensemble of transitions, which are saved to a file or returned as an array of paths. Multi-threading is enabled.

## Keyword arguments

  - `rad_i=0.1`: radius of ball around `x_i`
  - `rad_f=0.1`: radius of ball around `x_f`
  - `cut_start=true`: if `false`, returns the whole trajectory up to the transition
  - `Nmax`: number of attempts before the algorithm stops even if less than `N` transitions occurred.
  - `dt=0.01`: time step of integration
  - `tmax=1e3`: maximum time when the simulation stops even `x_f` has not been reached
  - `rad_dims=1:length(sys.u)`: the directions in phase space to consider when calculating the radii
    `rad_i` and `rad_f`. Defaults to all directions. To consider only a subspace of state space,
    insert a vector of indices of the dimensions to be included.
  - `solver=EM()`: numerical solver. Defaults to Euler-Mayurama
  - `progress`: shows a progress bar with respect to `Nmax`
  - `savefile`: if `nothing`, no data is saved to a file. To save to a file, see below.

See also [`transition`](@ref).

## Saving data to file

The `savefile` keyword argument allows saving the data to a `.jld2` or `.h5` file. To do so:

 1. Create and open a file by typing `file = jld2open("filename.jld2", "a+")` or `file = h5open("filename.h5", "cw")`. This requires `JLD2.jl`/`HDF5.jl`; the convenience functions [`make_jld2`](@ref), [`make_h5`](@ref) provide this out of the box.
 2. Pass the label `file` to the `savefile` argument of `transitions`.
 3. Don't forget to `close(file)` at the end.

## Output

`[samples, times, idx, N_fail]`

  - `samples` (Array of Matrices): sample paths. Each path i has size (dim × N_i), where N_i is the number of path points
  - `times` (Array of Vectors): time values (since simulation start) of path points for each path
  - `idx` (Array): list of sample indices i that produced a transition
  - `N_fail` (Int): number of samples that failed to produce a transition

> An example script using `transitions` is available [here](https://github.com/juliadynamics/CriticalTransitions.jl/blob/main/scripts/sample_transitions_h5.jl).
"""
function transitions(sys::CoupledSDEs, x_i, x_f, N = 1;
        rad_i = 0.1,
        rad_f = 0.1,
        tmax = 1e3,
        Nmax = 1000,
        cut_start = true,
        rad_dims = 1:length(current_state(sys)),
        savefile = nothing,
        showprogress::Bool = false,
        kwargs...)
    """
    Generates N transition samples of sys from x_i to x_f.
    Supports multi-threading.
    rad_i:      ball radius around x_i
    rad_f:      ball radius around x_f
    cut_start:  if false, saves the whole trajectory up to the transition
    savefile:   if not nothing, saves data to a specified open .jld2 file
    """

    samples, times, idx::Vector{Int64}, r_idx::Vector{Int64} = [], [], [], []

    # iterator = showprogress ? tqdm(1:Nmax) : 1:Nmax

    Threads.@threads for j in 1:Nmax
        sim, simt, success = transition(sys, x_i, x_f;
            rad_i = rad_i, rad_f = rad_f, rad_dims = rad_dims,
            tmax = tmax, cut_start = cut_start, kwargs...)

        if success
            if showprogress
                print("\rStatus: $(length(idx)+1)/$(N) transitions complete.")
            end

            if savefile == nothing
                push!(samples, sim)
                push!(times, simt)
            else # store or save in .jld2/.h5 file
                write(savefile, "paths/path " * string(j), sim)
                write(savefile, "times/times " * string(j), simt)
            end

            push!(idx, j)

            if length(idx) > max(1, N - Threads.nthreads())
                break
            else
                continue
            end
        else
            push!(r_idx, j)
        end
    end


    success_rate = length(idx) / (length(r_idx) + length(idx))
    mean_res_time = sum([times[i][1] for i in 1:length(times)]) + tmax * length(r_idx)
    mean_trans_time = mean([(times[i][end] - times[i][1]) for i in 1:length(times)])

    return TransitionPathEnsemble(
        samples, times, success_rate, mean_res_time, mean_trans_time)
end;
