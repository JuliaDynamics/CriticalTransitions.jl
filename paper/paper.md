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

Metastability and tipping phenomena are important features of nonlinear dynamical systems in the natural and human world. A key element of metastable behavior are critical transitions between distinct dynamical regimes, often triggered by random perturbations or time-dependent external forcing. Alongside an expanding mathematical theory of complex nonlinear dynamics, Julia software has become a leading code base to study dynamical systems numerically. However, functionality for analyzing critical behavior in stochastic and non-autonomous systems has been missing. CriticalTransitions.jl offers an intuitive, extendible and tested toolbox for critical transitions in forced dynamical systems. Built on a familiar interface for defining systems with stochastic and/or time-dependent forcing, the package allows users to sample, quantify and understand noise- and rate-induced transitions. Currently available features include: action minimization for computing most probable transition paths via large deviation theory, methods from transition path theory, and R-tipping concepts such as critical forcing rates.

# Statement of need

Critical transitions are a topic of growing scientific interest, given their relevance in diverse disciplines and applications. 

Solvers for nonautonomous and stochastic differential equations exist, but they have not been incorporated in the DynamicalSystems.jl interface.

# Concept

CriticalTransitions.jl addresses dynamical systems of the general form
$$ \text{d}\mathbf{x} = \mathbf{f}(\mathbf{x},\, p(t)) \,\text{d}t + \mathbf{g}(\mathbf{x},\, p(t))\, \text{d}\mathbf{\mathcal{N}}_t \,,$$
where the state $\mathbf{x}(t) \in \mathbb{R}^D$ evolves under the deterministic drift $\mathbf{f}$ and stochastic forcing described by a noise function $\mathbf{g}$ and noise process $\mathbf{\mathcal{N}}_t$. Both the drift and noise functions may depend explicitly on time via parameters $p$. This setup corresponds to the `CoupledSDEs` system type in DynamicalSystems.jl, based on DifferentialEquation.jl's `SDEProblem`. In the absence of noise, this system reduces to the `CoupledODEs` type (based on `ODEProblem`).

A widely studied case is  the noise function $\mathbf{g(\mathbf{x})} = \sigma \mathbf\Sigma (\mathbf x)$, where

The deterministic autonomous system $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}, p)$ 

Some problems are more conveniently described in a Hamiltonian formulation using generalized position and momentum coordinates. 

# bla

and refer to \autoref{eq:fourier} from text.

# Example

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

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

# Acknowledgements

We would like to thank Tobias Grafke, Oliver Mehling, Calvin Nesbitt and Jeroen Wouters for their input that helped develop CriticalTransitions.jl.

# References