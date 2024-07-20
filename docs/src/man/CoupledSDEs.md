# Define a CoupledSDE

A `CoupledSDEs` defines a stochastic dynamical system of the form

```math
\text{d}\vec x = f(\vec x(t); \ p)  \text{d}t + g(\vec x(t);  \ p) \text{d}\mathcal{W} \ ,
```
where $\text{d}\mathcal{W}=\Gamma \cdot \text{d}\mathcal{N}$, $\vec x \in \mathbb{R}^\text{dim}$ and $\mathcal N$ denotes a stochastic process. The (positive definite) noise covariance matrix is $\Sigma = \Gamma \Gamma^\top \in \mathbb R^{N\times N}$.

The function $f$ is the deterministic part of the system and is assumed to be of similar form as is accepted in [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/latest/tutorial/), i.e., `f(u, p, t)` for out-of-place (oop) and `f(du, u, p, t)` for in-place (iip).

The function $g$ represent the stochastics dynamics of the system and should be the of the same type (iip or oop) as $f$. The keyword `noise` defines the system [noise process](#noise-process). In combination with `g` one can define different type of stochastic systems. Examples of different type of stochastics systems can be found on the [StochasticDiffEq.jl tutorial page](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/). A quick overview of the different types of stochastic systems can be found [here](#Type-of-stochastic-system).

!!! info
    Note that nonlinear mixings of the Noise Process $\mathcal{W}$ are not Stochasitic Differential Equations but are a different class of differential equations of random ordinary differential equations (RODEs) which have a separate set of solvers. See [this example](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/rode_example/) of DifferentialEquations.jl.

!!! warning
    Some algrorithms need to know the overall noise strength of the system. However, the noise strength is ambiguous in the case of non-diagonal noise. Hence, these methods ask for the noise strength to be provided explicitly. This can be done by setting the `noise_strength` keyword.


```@docs
CoupledSDEs
```
## Type of stochastic system
Let us make some examples of the different types of stochastic systems that can be defined.
```@example type
using CriticalTransitions, Plots
import Random # hide
Random.seed!(10) # hide
f!(du, u, p, t) = du .= 1.01u
σ = 0.25
```
### Additive noise
When `g` is independent of the state variables `u`, the noise is called additive.

#### Diagonal noise
A system of diagional noise is the most common type of noise. It is defined by a vector of random numbers `dW` whose size matches the output of `g` where the noise is applied element-wise, i.e. `g.*dW`.
```@example type
t0 = 0.0; W0 = zeros(2);
W = WienerProcess(t0, W0, 0.0)
sde = CoupledSDEs(f!, diagonal_noise!(σ), zeros(2); noise=W)
```
or equivalently
```@example type
sde = CoupledSDEs(f!, diagonal_noise!(σ), zeros(2))
```
The `diagonal_noise!` function is a helper function equivalent to `(du, u, p, t) -> σ .* idfunc!(du, u, p, t)` with `idfunc!(du, u, p, t) = (du .= ones(length(u)); return nothing)`. The vector `dW` is by default zero mean white gaussian noise $\mathcal{N}(0, \text{d}t)$ where the variance is the timestep $\text{d}t$ unit variance (Wiener Process).
```@example type
sol = simulate(sde, 1.0, dt=0.01, alg=SOSRA())
plot(sol)
```

#### Scalar noise
Scalar noise is where a single random variable is applied to all dependent variables. To do this, one has to give the noise process to the `noise` keyword of the `CoupledSDEs` constructor. A common example is the Wiener process starting at `W0=0.0` at time `t0=0.0`.

```@example type
t0 = 0.0; W0 = 0.0;
noise = WienerProcess(t0, W0, 0.0)
sde = CoupledSDEs(f!, diagonal_noise!(σ), rand(2)./10; noise=noise)
sol = simulate(sde, 1.0, dt=0.01, alg=SOSRA())
plot(sol)
```
### Multiplicative noise
Multiplicative Noise is when $g_i(t, u)=a_i u$
```@example type
function g(du, u, p, t)
    du[1] = σ*u[1]
    du[2] = σ*u[2]
    return nothing
end
sde = CoupledSDEs(f!, g, rand(2)./10)
sol = simulate(sde, 1.0, dt=0.01, alg=SOSRI())
plot(sol)
```

#### Non-diagonal noise
Non-diagonal noise allows for the terms to linearly mixed via g being a matrix. Suppose we have two Wiener processes and two dependent random variables such that the output of `g` is a 2x2 matrix. Therefore, we have
```math
du_1 = f_1(u,p,t)dt + g_{11}(u,p,t)dW_1 + g_{12}(u,p,t)dW_2 \\
du_2 = f_2(u,p,t)dt + g_{21}(u,p,t)dW_1 + g_{22}(u,p,t)dW_2
```
To indicate the structure that `g` should have, we can use the `noise_rate_prototype` keyword. Let us define a special type of non-diagonal noise called commutative noise. For this we can utilize the `RKMilCommute` algorithm which is designed to utilise the structure of commutative noise.

```@example type
function g(du, u, p, t)
  du[1,1] = σ*u[1]
  du[2,1] = σ*u[2]
  du[1,2] = σ*u[1]
  du[2,2] = σ*u[2]
    return nothing
end
sde = CoupledSDEs(f!, g, rand(2)./10, noise_rate_prototype = zeros(2, 2))
sol = simulate(sde, 1.0, dt=0.01, alg=RKMilCommute())
plot(sol)
```

!!! warning
    Non-diagonal problem need specific type of solvers. See the [SciML recommendations](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#sde_solve).

### Corelated noise
```@example type
ρ = 0.3
Σ = [1 ρ; ρ 1]
t0 = 0.0; W0 = zeros(2); Z0 = zeros(2);
W = CorrelatedWienerProcess(Σ, t0, W0, Z0)
sde = CoupledSDEs(f!, diagonal_noise!(σ), zeros(2); noise=W)
sol = simulate(sde, 1.0, dt=0.01, alg=SOSRA())
plot(sol)
```

## Noise process
We provide the noise processes $\text{d}\mathcal{W}$ that can be used in the stochastic simulations through the [DiffEqNoiseProcess.jl](https://docs.sciml.ai/DiffEqNoiseProcess/stable) package. A complete list of the available processes can be found [here](https://docs.sciml.ai/DiffEqNoiseProcess/stable/noise_processes/). We list some of the most common ones below:
```@docs
WienerProcess
SimpleWienerProcess
OrnsteinUhlenbeckProcess
CorrelatedWienerProcess
```