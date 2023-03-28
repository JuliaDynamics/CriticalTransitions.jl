"""
    equilib(sys::StochSystem, state::State; kwargs...)
Returns the equilibrium solution of the system `sys` for given initial condition `state`.

> Warning: This algorithm simply evolves the deterministic system forward in time until a steady-state condition is satisfied.
> Thus, the algorithm may output a false solution if it gets stuck in a quasi-equilibrium, or slowly evolving state.
> For more robust results, use `fixedpoints`.

## Keyword arguments:
* `abstol = 1e-5`: steady-state condition. Simulation ends when the rate of change (Euclidean distance in state space) of the state falls below `abstol`.
* `tmax = 1e5`: maximum simulation time before the algorithm stops even if the steady-state condition is not reached.
* `dt = 0.01`: time step of the ODE solver.
* `solver = Euler()`: ODE solver used for evolving the state.
"""
function equilib(sys::StochSystem, state::State;
    dt=0.01,
    tmax=1e5,
    abstol=1e-5,
    solver=Euler())
    
    condition(u, t, integrator) = norm(integrator.uprev-u) < abstol
    affect!(integrator) = terminate!(integrator)
    equilib_cond = DiscreteCallback(condition, affect!)

    prob = ODEProblem(sys.f, state, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=equilib_cond, save_on=false, save_start=false)

    sol.u[1]
end;

function fixedpoints(sys::StochSystem, box)
    jac(u,p,t) = ForwardDiff.jacobian((x) -> sys.f(x,p,t), u)
    DynamicalSystems.fixedpoints(CoupledODEs(sys), box, jac)
end;

"""
    fixedpoints(sys::StochSystem, bmin::Vector, bmax::Vector)
Returns fixed points, their eigenvalues and stability of the system `sys` within the state space volume defined by `bmin` and `bmax`.

> This is a wrapper around the [`fixedpoints`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/chaos/periodicity/#ChaosTools.fixedpoints) function of `DynamicalSystems.jl`.

## Input
* `bmin` (Vector): lower limits of the state space box to be considered, as a vector of coordinates
* `bmax` (Vector): upper limits
* alternatively `box` (IntervalBox) can replace `bmin` and `bmax`

> Example: `fixedpoints(sys, [-2,-1,0], [2,1,1])` finds the fixed points of the 3D system `sys` in a cube defined by the intervals `[-2,2] × [-1,1] × [0,1]`.

## Output
`[fp, eigs, stable]`
* `fp`: `Dataset` of fixed points
* `eigs`: vector of Jacobian eigenvalues of each fixed point
* `stable`: vector of booleans indicating the stability of each fixed point (`true`=stable, `false`=unstable)

## Additional methods
* `fixedpoints(sys::StochSystem, box)`
"""
function fixedpoints(sys::StochSystem, bmin::Vector, bmax::Vector)
    box = intervals_to_box(bmin, bmax)
    jac(u,p,t) = ForwardDiff.jacobian((x) -> sys.f(x,p,t), u)
    DynamicalSystems.fixedpoints(CoupledODEs(sys), box, jac)
end;

function saddles_idx(fps::Tuple)
    num = size(fps[1],1); # number of fixed points
    dim = size(fps[1],2); # dimension of the system
    eigenvalues = fps[2];
    idx = [false for i ∈ 1:num];
    for ii ∈ 1:num
        imag_parts = [imag(eigenvalues[ii][jj]) for jj ∈ 1:dim]
        if all(imag_parts.==0) # we have purely real eigenvalues
            real_parts = [real(eigenvalues[ii][jj]) for jj ∈ 1:dim];
            if prod(real_parts) < 0 # we have at least positive eigenvalue and at least one negative eigenvalue
                idx[ii] = true;
            end
        end
    end
    idx
end

function repellers_idx(fps::Tuple)
    num = size(fps[1],1); # number of fixed points
    dim = size(fps[1],2); # dimension of the system
    eigenvalues = fps[2];
    idx = [false for i ∈ 1:num];
    for ii ∈ 1:num
        real_parts = [real(eigenvalues[ii][jj]) for jj ∈ 1:dim];
        if all(real_parts .> 0) 
            idx[ii] = true;
        end
    end
    idx
end

function attractors_idx(fps::Tuple)
    num = size(fps[1],1); # number of fixed points
    dim = size(fps[1],2); # dimension of the system
    eigenvalues = fps[2];
    idx = [false for i ∈ 1:num];
    for ii ∈ 1:num
        real_parts = [real(eigenvalues[ii][jj]) for jj ∈ 1:dim];
        if all(real_parts .< 0) 
            idx[ii] = true;
        end
    end
    idx
end