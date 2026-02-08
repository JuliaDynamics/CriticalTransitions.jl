---
title: 'CriticalTransitions.jl: A toolbox for noise- and rate-induced transitions in dynamical systems'
tags:
  - Julia
  - dynamical systems
  - critical transitions
  - stochastic dynamics
  - large deviation theory
authors:
  - name: Reyk Börner
    orcid: 0000-0003-4593-7721
    equal-contrib: true
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
affiliations:
 - name: Institute for Marine and Atmospheric research Utrecht, Utrecht University, The Netherlands
   index: 1
date: 1 July 2025
bibliography: paper.bib
---

# Summary

Metastability and tipping phenomena are important features of nonlinear dynamical systems in the natural and human world. A key element of metastable behavior are critical transitions between distinct dynamical regimes, often triggered by random perturbations or time-dependent external forcing. The Julia package CriticalTransitions.jl offers a numerical toolbox for simulating and studying critical behavior, particularly noise- and rate-induced transitions. Alongside recent theoretical and methodological advances for forced dynamical systems, Julia software has become a leading code base for fast, reliable and user-friendly implementations of dynamical systems theory. This package builds on DynamicalSystems.jl and DifferentialEquations.jl to add intuitive, extendable and tested functionality for stochastic and non-autonomous systems. Focused on sampling, quantifying and understanding transitions, available features include: action minimization for computing most probable transition paths via large deviation theory, methods from transition path theory, and R-tipping concepts such as critical forcing rates.

# Statement of need

Critical transitions are a topic of growing scientific interest, given their relevance in diverse disciplines and applications -- from protein folding to climate tipping points. Research problems are often not tractable analytically, requiring numerical techniques for applying ideas from dynamical systems theory. Since transition events are typically rare, efficient numerics are essential.

By combining speed and user-friendliness, Julia is a promising language for making these tools accessible. While SciML.jl and DifferentialEquations.jl provide fast, powerful solvers for nonautonomous and stochastic systems, DynamicalSystems.jl offers a user interface to construct the systems and study their autonomous, deterministic dynamics. However, a similar interface and functionality for forced dynamical systems has been missing. Specifically, various powerful algorithms have been proposed to sample and predict transition paths but, to our knowledge, had not been implemented in Julia (or, sometimes, in any common programming language) in a generic way. CriticalTransitions.jl fills this gap and is designed to grow with future contributions.

# Concept

CriticalTransitions.jl addresses dynamical systems of the general form
$$ \text{d}\bm{x} = \bm{f}(\bm{x},\, p(t)) \,\text{d}t + \bm{g}(\bm{x},\, p(t))\, \text{d}\bm{\mathcal{N}}_t \,,$$
where the state $\bm{x}(t) \in \mathbb{R}^D$ evolves under the deterministic drift $\bm{f}$ and stochastic forcing described by a noise function $\bm{g}$ and noise process $\bm{\mathcal{N}}_t$. Both the drift and noise functions may depend explicitly on time via parameters $p$. This setup corresponds to the `CoupledSDEs` system type in DynamicalSystems.jl, based on DifferentialEquation.jl's `SDEProblem`. In the absence of noise, this system reduces to the `CoupledODEs` type (based on `ODEProblem`).

Much research is devoted to dynamical systems driven by weak Gaussian noise, in which transitions and other rare events can be studied via large deviation theory. Here the noise function becomes $\bm{g(\bm{x})} = \sigma \bm\Sigma (\bm x)$, where the noise strength $\sigma$ and diffusion matrix $\bm \Sigma$ can be conveniently defined as keyword arguments of the `CoupledSDEs` constructor. For this class of systems, CriticalTransitions.jl provides functionality to calculate action functionals and most probable transition paths (instantons) by numerically solving the associated action minimization problem. Some problems are more conveniently described in a Hamiltonian formulation using generalized position and momentum coordinates, for which the `ExtendedHamilton` type is provided. Langevin dynamics can also be formulated as a `LangevinSystem`.

Another broadly studied class are nonautonomous dynamical systems where the (deterministic) dynamics are driven by a time-dependent modulation of the system parameters: $\dot{\bm{x}} = \bm{f}(\bm{x}, p(t))$. Here, the question is whether the forcing protocol $p(t)$ induces critical behavior in the responding system. In CriticalTransitions.jl, the forcing can conveniently be specified as a `RateConfig` type and applied to the system of interest via `apply_ramping`. Combined with the functionality of the Attractors.jl package, this interface forms the basis to study rate-induced tipping (R-tipping).

In many real-world applications, stochastic and parametric forcings act in parallel, and transitions are often caused by the combined effect of random fluctuations and external driving. By connecting both aspects in CriticalTransitions.jl, the package enables extensions beyond the classical problem settings of large deviation theory on the one hand and systems subjected to parameter drift on the other.

 
# Example
To illustrate some key functionality, let us consider the bistable FitzHugh-Nagumo model, a two-dimensional system of ordinary differential equations originating from theoretical neuroscience (see @borner2024saddle):
$$\begin{align}
\frac{\text d x}{\text dt} &= \nu\big[- x^3 + x - y + I(t) \big] \nonumber \\
\frac{\text d y}{\text dt} &= - \beta y + x \,.
\end{align}$$
The evolution of the system state $\bm x(t) \equiv (x(t), y(t))$ depends on the parameters $p=(\nu, \beta, I)$, where $\nu=0.1$ sets the timescale separation between $x$ and $y$, and $I$ denotes an external driving. For $I=0$ and $\beta=3$, the system possesses two stable fixed points at $\bm x_{A, B} = \pm (xx, yy)$. After discussing the stochastically forced case with $I=0$, we also consider non-autonomous forcing $I(t)$.

First, suppose the system is driven by Gaussian white noise with noise strength $\sigma=0.1$ and covariance matrix $\bm \Sigma\bm \Sigma^\top = \begin{pmatrix}1 & 1/3 \\ 1/3 & 2\end{pmatrix}$. We can easily set this problem up as a `CoupledSDEs`:
```julia
using CriticalTransitions

function fitzhugh_nagumo(u, p, t)
  x, y = u
  nu, beta, I = p

  dx = nu*(-x^3 + x - y + I)
  dy = -beta*y + x

  return SVector{2}([dx, dy])
end

p = [0.1, 3, 0]
initial_condition = [1.0, 0.0]
noise_strength = 0.1
covariance = [1 1/3; 1/3 2]

sys = CoupledSDEs(fitzhugh_nagumo, initial_condition, p; noise_strength, covariance)
```

Starting from $\bm x_A$, how long does it typically take until trajectories transition to the competing attractor $\bm x_B$? What path do sampled noise-induced transitions take, and what is the most probable transition path in the weak-noise limit ($\sigma\to 0$)?

Ensembles of $N$ transition paths from $\bm x_A$ to (a neighborhood around) $\bm x_B$ can be generated via parallelized Monte Carlo rejection sampling with the `transitions` function, which also returns transition statistics such as their rareness.

```julia
x_A = []; x_B = []
N = 100 # sample size

transitionsAB = transitions(sys, x_A, x_B, N)
```
CriticalTransitions.jl provides several methods to compute most probable transition paths (instantons). For example, the geometric minimum action method (which minimizes the Freidlin-Wentzell action in a time-independent formulation) is accessible as follows:

```julia
instantonAB = geometric_min_action(sys, x_A, x_B; abstol=1e-4)
```
Here we instructed the algorithm to stop once the action value changes by less than $10^{-4}$ between iterations of the minimization procedure. The resulting minimum action path and ensemble of sample transition paths is shown in Fig. 1.

Now, suppose there is no noise ($\sigma=0$) but the external forcing parameter $I$ varies in time, changing nonlinearly from $I_\text{min}=-0.185$ to $I_\text{max}=0.185$ according to
$$\begin{align} I(t) = \frac{I_\text{max}-I_\text{min}}{2} \big(\tanh(rt) + 1\big) + I_\text{min} \,, \end{align}$$
where the rate parameter $r$ determines how rapidly the forcing changes relative to the internal timescales of the FitzHugh-Nagumo model. We must thus turn the stochastic FitzHugh-Nagumo model defined as `sys` above into a non-autonomous, deterministic version `sys_ramp`. We can easily achieve this by first converting to a deterministic system (`CoupledODEs`) and then applying the forcing protocol (`RateConfig`) to make the parameter $I$ time-dependent.

```julia
sys_deterministic = CoupledODEs(sys)

I_min, I_max = -0.185, 0.185
p_ramping = [8.0] # rate parameter

I(t, p_ramping) = (I_max - I_min)/2 * (tanh(p_ramping[1]*t) + 1) + I_max

ramp = RateConfig(...)
sys_ramp = apply_ramping(sys_deterministic, ramp)
```

One important question is: how does the time-dependent forcing change the equilibria of the system? With CriticalTransitions.jl, this can be directly analyzed using the `moving_sinks` function, which return the result shown in Fig. 2.

```julia
equilib = moving_sinks(...)
```

Within the given forcing range, the model remains bistable but the positions of the stable fixed points change as a function of $I(t)$. If we initialize the system at $\bm x_A$ at time $t=-100$, will the system keep tracking that fixed point or transition to $\bm x_B$ as a result of the parameter change? This depends on the forcing rate $r$: here, if $r<r_c$, the trajectory tracks $\bm x_A$, whereas for $r>r_c$ the system undergoes rate-induced tipping. We can compute the critical rate $r_c$ via the `critical_rate` function:

```julia
```
The results are displayed in Fig. 3.

Lastly, we can also use CriticalTransitions.jl to study the combined effect of noise and parameter ramping. For example, if we apply the Gaussian stochastic forcing from before to the ramped system near its critical rate, we observe that sample transition paths exhibit ensemble splitting --  where some ensemble members transition and others do not within a given time interval.

# Further reading
The full functionality of CriticalTransitions.jl is documented at [juliadynamics.github.io/CriticalTransitions.jl](https://juliadynamics.github.io/CriticalTransitions.jl/stable/), including various code examples.

# Acknowledgements

We would like to thank Tobias Grafke, Oliver Mehling, Calvin Nesbitt and Jeroen Wouters for their input that helped develop CriticalTransitions.jl.

# References

# Instructions

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

and refer to \autoref{eq:fourier} from text.


For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Another section

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
