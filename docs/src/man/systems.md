# Predefined dynamical systems

## FitzHugh-Nagumo system

A 2D system given by

```math
\begin{aligned}
\frac{du}{dt} &= \frac{1}{\epsilon} \left( -\alpha u^3 + \gamma u - \kappa v + I \right) \\
\frac{dv}{dt} &= -\beta v + u \ ,
\end{aligned}
```

where ``\epsilon`` is the parameter of time scale separation between the state variables ``u`` and ``v``. The parameters ``\alpha >0``, ``\beta >1``, ``\gamma>0``, and ``\kappa>0`` are real constants, and ``I`` denotes a driving term.

```@docs
FitzHughNagumo(u, p, t)
FitzHughNagumo!(du, u, p, t)
```

## Modified Truscott-Brindley system

A two-dimensional predator-prey system given by 

```math
\begin{aligned}
\frac{\text{d} P}{\text{d} t} &= rP\bigg(1-\frac{P}{K}\bigg)-a\frac{P^2}{h^2+P^2}Z \\
\frac{\text{d} Z}{\text{d} t} &= \xi\bigg(\frac{aP^2}{h^2+P^2}Z - mZ^2\bigg)
\end{aligned}
```
The variables ``P`` and ``Z`` represent the concentration of phytoplankton (prey) and zooplankton (predator) species respectively. The system parameters can be interpreted as follows: 

* ```r```: the phytoplankton growth rate
* ```K```: the carrying capacity of the phytoplankton
* ```a```: the maximal predation rate
* ```h```: the half-saturation constant
* ```\xi```: the time-scale separation between the evolution of phytoplankton and zooplankton 
* ```m```: the zooplankton mortality rate 

The model (cite...) is a modified version of the Truscott-Brindley system (cite ...), with a quadratic mortality term to enable bistability. In non-dimensional form, the above dynamics can be written...       

## Thermohaline Circulation box models

```@docs
rooth_smooth(u, p, t; kwargs...)
```
