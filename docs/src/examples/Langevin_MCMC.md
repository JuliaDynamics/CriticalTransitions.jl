```@meta
EditURL = "../../../examples/Langevin_MCMC.jl"
```

# Langevin Markov Chain Monte Carlo (LMCMC)

The LMCMC

## Required Packages

````@example LMCMC
````

## Problem Setup

We'll consider a 2D system with a double-well flow, called the Maier-Stein model. It is a famous benchmark problem as it exhibits non-gradient dynamics with two stable equilibrium points at (-1,0) and (1,0), connected by a non-trivial transition path.
The system's deterministic dynamics are given by:

Define the vector field

````@example LMCMC
f(u, v) = [u - u^3 - 10*u*v^2,  -(1 - u^2)*v]
f(x) = f(x...)
nothing # hide
````

## Optimal Control Formulation

The minimal action path minimizes the deviation from the deterministic dynamics:

````@example LMCMC
function ocp(T)
    action = @def begin
        t ∈ [0, T], time
        x ∈ R², state
        u ∈ R², control
        x(0) == [-1, 0]    # Starting point (left well)
        x(T) == [1, 0]     # End point (right well)
        ẋ(t) == u(t)       # Path dynamics
        ∫( sum((u(t) - f(x(t))).^2) ) → min  # Minimize deviation from deterministic flow
    end
    return action
end
nothing # hide
````
