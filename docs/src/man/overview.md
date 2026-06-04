# Choosing an approach

CriticalTransitions.jl provides several complementary methods that all answer the same underlying question: **how, where, and how often does a system undergo a rare transition between metastable states?** They differ in the regime where they are valid (noise strength, state-space dimension) and in what they return (a scalar rate, the most probable path, a spatial flux field). This page is the map; each method has its own manual page with the details.

## The method menu

| Your goal | Reach for | Regime / output |
|---|---|---|
| Observe transitions directly | [Direct sampling](@ref "Sampling transitions") | any noise level; brute-force trajectories |
| Most probable path **and** exponential rate | [Large deviation theory](@ref "Large deviation theory") | weak noise (``\varepsilon \to 0``); instanton & quasipotential |
| Escape rate, mean first-passage time, (quasi-)stationary density | [Rates, distributions & the generator](@ref "Rates, distributions & the generator") | finite noise; spectral observables |
| Where reactive paths concentrate, transition channels, committor | [Transition path theory](@ref "Transition Path Theory") | finite noise; spatial flux |
| Transition driven by parameter drift rather than noise | [Rate-induced transitions](@ref "Rate-induced transitions") | non-autonomous forcing |

## Picking by what you need

- **Only a rate or waiting time?** Use the spectral observables (`mean_first_passage_time`, quasi-stationary distribution) on the [generator](@ref "Rates, distributions & the generator"), or Kramers asymptotics for high barriers.
- **The transition path itself?** Use [large deviation theory](@ref "Large deviation theory") (the instanton via MAM / gMAM / sgMAM / shooting).
- **Spatial flux and channels?** Use [transition path theory](@ref "Transition Path Theory").
- **High barrier, weak noise?** Direct sampling becomes exponentially expensive; prefer LDT or the spectral observables.
- **Shallow barrier?** Just simulate with [direct sampling](@ref "Sampling transitions").
- **Parameter drift, not noise?** Use [rate-induced transitions](@ref "Rate-induced transitions").

## Composing methods

The methods are not mutually exclusive. Transition path theory and the spectral rate observables are both built on the same discretised [generator](@ref "Rates, distributions & the generator"). In large deviation theory, a fast gMAM / sgMAM solve makes a good warm start for the more delicate multiple-shooting boundary-value problem; see [Large deviation theory](@ref "Large deviation theory").
