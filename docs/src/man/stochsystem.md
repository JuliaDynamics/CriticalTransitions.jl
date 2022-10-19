# Defining a StochSystem

A `StochSystem` defines a stochastic dynamical system of the form

``\text{d}\vec x = f(\vec x(t); \ p_f)  \text{d}t + \sigma g(\vec x(t);  \ p_g)  \Sigma \cdot \text{d}\mathcal{N} \ ,``

where $\vec x \in \mathbb{R}^\text{dim}$ and $\mathcal N$ denotes a stochastic process.

An instance of StochSystem is created via `StochSystem(f, pf, dim [, σ [, g, pg, Σ , process]])`,
taking the following arguments:
* `f` (Function): Dynamical ODE rule describing the drift term of the system, corresponding to `f` in the ODEProblem of `DifferentialEquations`. Can be defined in-place (`f!(du, u, p, t)`) or out-of-place (`f(u,p,t)`).
* `pf` (Vector or Nothing): Parameter vector for the drift term.
* `dim` (Int64): Dimension of the system.
* `σ` (Float64): Noise intensity. Defaults to `1.0`.
* `g` (Function): Dynamical ODE rule describing the noise term of the system. Same format as `f`. Defaults to `idfunc`.
* `pg` (Vector or Nothing): Parameter vector of the noise term.
* `Σ` (Matrix): Noise covariance matrix. Defaults to `I` (identity matrix).
* `process` (String): Noise process. Defaults to `white-gauss` (independent n-dimensional Brownian motion).

### Shortcut methods
* `StochSystem(f, pf, dim)`
* `StochSystem(f, pf, dim, σ)`
* `StochSystem(f, pf, dim, σ, Σ)`

```@docs
StochSystem(f, pf, dim [, σ [, g, pg, Σ, process]])
```