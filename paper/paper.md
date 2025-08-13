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

Metastability and tipping phenomena are important features of nonlinear dynamical systems, occurring in various contexts in the natural and human world. A key element of metastable behavior are critical transitions between distinct dynamical regimes, often triggered by random perturbations or time-dependent external forcing. Alongside an expanding mathematical theory of complex nonlinear dynamics, Julia software has become a leading code base to study dynamical systems numerically. However, functionality for analyzing critical behavior in stochastic and non-autonomous systems has been missing. CriticalTransitions.jl offers a user-friendly, coherent, extendible and tested toolbox for critical transitions in forced dynamical systems. Built on an intuitive interface for defining systems with stochastic and/or time-dependent forcing, the package allows users to sample, quantify and understand noise- and rate-induced transitions. Available features include rare event sampling, computing most probable transition paths using large deviation theory, methods from transition path theory, and determining critical forcing rates. Both the source code and user interface are designed to fit seamlessly to DynamicalSystems.jl.

# Statement of need


# Concept

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
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

We acknowledge contributions from Tobias Grafke, Oliver Mehling, Calvin Nesbitt and Jeroen Wouters.

# References