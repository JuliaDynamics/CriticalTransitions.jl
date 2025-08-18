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

Metastability and tipping phenomena are important features of nonlinear dynamical systems in the natural and human world. A key element of metastable behavior are critical transitions between distinct dynamical regimes, often triggered by random perturbations or time-dependent external forcing. The Julia package CriticalTransitions.jl offers a numerical toolbox for simulating and studying critical behavior, particularly noise- and rate-induced transitions. Alongside recent theoretical and methodological advances for forced dynamical systems, Julia software has become a leading code base for fast, reliable and user-friendly implementations of dynamical systems theory. This package builds on DynamicalSystems.jl and DifferentialEquations.jl to add intuitive, extendible and tested functionality for stochastic and non-autonomous systems. Focused on sampling, quantifying and understanding transitions, available features include: action minimization for computing most probable transition paths via large deviation theory, methods from transition path theory, and R-tipping concepts such as critical forcing rates.

# Statement of need

Critical transitions are a topic of growing scientific interest, given their relevance in diverse disciplines and applications -- from protein folding to climate tipping points. Research problems are often not tractable analytically, requiring numerical techniques for applying ideas from dynamical systems theory. Since transition events are typically rare, efficient numerics are essential.

By combining speed and user-friendliness, Julia is a promising language for making these tools accessible. While SciML.jl and DifferentialEquations.jl provide fast, powerful solvers for nonautonomous and stochastic systems, DynamicalSystems.jl offers a user interface to construct the systems and study their autonomous, deterministic dynamics. However, a similar interface and functionality for forced dynamical systems had been missing. Specifically, various powerful algorithms have been proposed to sample and predict transition paths but, to our knowledge, had not been implemented in Julia (or, sometimes, in any common programming language) in a generic way. CriticalTransitions.jl fills this gap and is designed to grow with future contributions.

# Concept

CriticalTransitions.jl addresses dynamical systems of the general form
$$ \text{d}\mathbf{x} = \mathbf{f}(\mathbf{x},\, p(t)) \,\text{d}t + \mathbf{g}(\mathbf{x},\, p(t))\, \text{d}\mathbf{\mathcal{N}}_t \,,$$
where the state $\mathbf{x}(t) \in \mathbb{R}^D$ evolves under the deterministic drift $\mathbf{f}$ and stochastic forcing described by a noise function $\mathbf{g}$ and noise process $\mathbf{\mathcal{N}}_t$. Both the drift and noise functions may depend explicitly on time via parameters $p$. This setup corresponds to the `CoupledSDEs` system type in DynamicalSystems.jl, based on DifferentialEquation.jl's `SDEProblem`. In the absence of noise, this system reduces to the `CoupledODEs` type (based on `ODEProblem`).

Much research is devoted to dynamical systems driven by weak Gaussian noise, in which transitions and other rare events can be studied via large deviation theory. Here the noise function becomes $\mathbf{g(\mathbf{x})} = \sigma \mathbf\Sigma (\mathbf x)$, where the noise strength $σ$ and diffusion matrix $\mathbf \Sigma$ can be conveniently defined as keyword arguments of the `CoupledSDEs` constructor. For this class of systems, CriticalTransitions.jl provides functionality to calculate action functionals and most probable transition paths (instantons) by numerically solving the associated action minimization problem. Some problems are more conveniently described in a Hamiltonian formulation using generalized position and momentum coordinates, for which the `ExtendedHamilton` type is provided. Langevin dynamics can also be formulated as a `LangevinSystem`.

Another broadly studied class are nonautonomous dynamical systems where the (deterministic) dynamics are driven by a time-dependent modulation of the system parameters: $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}, p(t))$. Here, the question is whether the forcing protocol $p(t)$ induces critical behavior in the responding system. In CriticalTransitions.jl, the forcing can conveniently be specified as a `RateConfig` type and applied to the system of interest via `apply_ramping`. Combined with the functionality of the Attractors.jl package, this interface forms the basis to study rate-induced tipping (R-tipping).

In many real-world applications, stochastic and parametric forcings act in parallel, and transitions are often caused by the combined effect of random fluctuations and external driving. By connecting both aspects in CriticalTransitions.jl, the package enables extensions beyond the classical problem settings of large deviation theory on one hand and systems subjected to parameter drift on the other.

 
# Example
To illustrate some of the key functionality, let us consider the bistable FitzHugh-Nagumo model, a conceptual 2-dimensional system originating from theoretical neuroscience (see @borner2024saddle):
$$\begin{align}
\dot u &= - u^3 + u - v + I \nonumber \\
\dot v &= - \beta v + u
\end{align}$$
For $\beta=3$, 

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
