# Define a CoupledSDE

A `CoupledSDEs` defines a stochastic dynamical system of the form

``\text{d}\vec x = f(\vec x(t); \ p)  \text{d}t + g(\vec x(t);  \ p) \text{d}\mathcal{W} \ ,``
where $\text{d}\mathcal{W}=\Gamma \cdot \text{d}\mathcal{N}$, $\vec x \in \mathbb{R}^\text{dim}$ and $\mathcal N$ denotes a stochastic process. The (positive definite) noise covariance matrix is $\Sigma = \Gamma \Gamma^\top \in \mathbb R^{N\times N}$.

The function $f$ is the deterministic part of the system and is assumed to be of similar form as is accepted in [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/latest/tutorial/), i.e., `f(u, p, t)` for out-of-place (oop) and `f(du, u, p, t)` for in-place (iip).

The function $g$ represent the stochastics dynamics of the system and should be the of the same type (iip or oop) as $f$. The keyword `noise` defines the system [noise process](#noise-process). In combination with `g` one can define different type of stochastic systems. Examples of different type of stochastics systems can be found on the [StochasticDiffEq.jl tutorial page](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/). A quick overview of the different types of stochastic systems can be found [here](#Type-of-stochastic-system).

!!! warning
    Some algrorithms need to know the overall noise strength of the system. Howeverm the noise strength is ambiguous in the case of non-diagonal noise. Hence, these methods ask for the noise strength to be provided explicitly. This can be done by setting the `noise_strength` keyword.


```@docs
CoupledSDEs
```
## Type of stochastic system

```@example type
using CriticalTransitions

function meier_stein!(du, u, p, t) # out-of-place
    x, y = u
    du[1] = x - x^3 - 10 * x * y^2
    dv[2] = -(1 + x^2) * y
    return
end
σ = 0.1
```

### Diagonal noise
that is a vector of random numbers dW whose size matches the output of g where the noise is applied element-wise,
### scalar noise
 scalar noise where a single random variable is applied to all dependent variables
### Non-diagonal noise
 more general type of noise allows for the terms to linearly mixed via g being a matrix.

 In our `g` we define the functions for computing the values of the matrix.
We can now think of the SDE that this solves as the system of equations

```math
du_1 = f_1(u,p,t)dt + g_{11}(u,p,t)dW_1 + g_{12}(u,p,t)dW_2 \\
du_2 = f_2(u,p,t)dt + g_{21}(u,p,t)dW_1 + g_{22}(u,p,t)dW_2
```


!!! info
    Note that nonlinear mixings are not Stochasitic Differential Equations but are a different class of differential equations of random ordinary differential equations (RODEs) which have a separate set of solvers. See this example of [DiffernetialEquations.jl](@exref rode_example).
### Corelated noise

## Noise process
We provide the noise processes $text{d}\mathcal{W}$ that can be used in the stochastic simulations through the [DiffEqNoiseProcess.jl](https://docs.sciml.ai/DiffEqNoiseProcess/stable) package. A complete list of the available processes can be found [here](https://docs.sciml.ai/DiffEqNoiseProcess/stable/noise_processes/). We list some of the most common ones below:
