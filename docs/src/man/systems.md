# Predefined dynamical systems

## Life sciences
### FitzHugh-Nagumo system

A 2D system given by

```math
\begin{aligned}
\frac{du}{dt} &= \frac{1}{\epsilon} \left( -\alpha u^3 + \gamma u - \kappa v + I \right) \\
\frac{dv}{dt} &= -\beta v + u \ ,
\end{aligned}
```

where ``\epsilon`` is the parameter of time scale separation between the state variables ``u`` and ``v``. The parameters ``\alpha >0``, ``\beta >1``, ``\gamma>0``, and ``\kappa>0`` are real constants, and ``I`` denotes a driving term.

```@docs
fitzhugh_nagumo(u, p, t)
fitzhugh_nagumo!(du, u, p, t)
fhn_ϵσ(ϵ, σ)   
```

## Population dynamics
### Modified Truscott-Brindley system

A two-dimensional predator-prey system given by 

```math
\begin{aligned}
\frac{\text{d} P}{\text{d} t} &= rP\bigg(1-\frac{P}{K}\bigg)-\frac{aP^2}{h^2+P^2}Z \\
\frac{\text{d} Z}{\text{d} t} &= \xi\bigg(\frac{aP^2}{h^2+P^2}Z - mZ^2\bigg)
\end{aligned}
```
The variables ``P`` and ``Z`` represent the concentration of phytoplankton (prey) and zooplankton (predator) species respectively. The system parameters can be interpreted as follows: 

* ``r = 1``: the phytoplankton growth rate
* ``K = 1``: the carrying capacity of the phytoplankton
* ``a = 1/9``: the maximal predation rate
* ``h = 5/112``: the half-saturation constant
* ``\xi = 1/10``: the time-scale separation between the evolution of phytoplankton and zooplankton 
* ``m = 0.0525``: the zooplankton mortality rate 

The model (as discussed [here](http://dx.doi.org/10.1016/j.ecocom.2014.10.003)) is a modified version of the original 1994 Truscott-Brindley system (see [here](https://doi.org/10.1007/BF02458277)), with a quadratic mortality term to enable bistability. The variable ``t`` has units of days. Following a non-dimensionalisation of the system, the dynamical equations convert to:

```math
\begin{aligned}
\displaystyle\frac{\text{d} \tilde{P}}{\text{d} \tau} &= \alpha \tilde{P}\bigg(1-\beta\frac{\tilde{P}}{P_1}\bigg)-\frac{\gamma (\tilde{P}/\sqrt{P_1})^2(\tilde{Z}/Z_1)}{1+(\tilde{P}/P_1)^2} \\
\displaystyle\frac{\text{d} \tilde{Z}}{\text{d} \tau} &= \xi\tilde{Z}\bigg(\frac{(\tilde{P}/P_1)^2}{1+(\tilde{P}/P_1)^2} - \frac{\tilde{Z}}{Z_1}\bigg)
\end{aligned}
```

Specifically, here we introduce the dimensionless variables ``\tilde{P}, \tilde{Z}`` and ``\tau`` according to the following transformations: 
```math
\begin{aligned}
\frac{\tilde{P}}{P_1} = \frac{P}{P_0}, && \frac{\tilde{Z}}{Z_1} = \frac{Z}{Z_0}, && \tau = \frac{t}{t_0}
\end{aligned}
```  
where ``P_1, P_0, Z_1, Z_0`` and ``t_0`` are constants to be determined in the non-dimensionalisation process. We take ``P_0 = h, Z_0 = a/m, t_0 = 1/a`` and subsequently introduce the following dimensionless parameters ``\alpha := r/a, \beta := h/K, \gamma := a/(mh)``. The values of the (dimensionless) constants ``P_1`` and ``Z_1`` are tuned according to the user's preference; they can be chosen in  a way such that the fixed points of the system are all contained within the ``[0,1]\times [0,1]`` subspace, for instance (for which the values ``P_1 = \beta`` and ``Z_0 = 5/6`` are appropriate in this set up).  

Below, the functions [`modifiedtruscottbrindley`](@ref) and [`modifiedtruscottbrindley!`](@ref) implement the non-dimensional form of the system. 

```@docs
modifiedtruscottbrindley(u, p, t)
modifiedtruscottbrindley!(du, u, p, t)
modtb_αξσ(α, ξ, σ)
```

## Earth & Climate
### Thermohaline Circulation box models

#### Stommel's hemispheric 2-box model
The famous 2-box model introduced by Stommel (1961) and later studied by Cessi (1994) can be generalized in the following form,

```math
\begin{aligned}
\frac{dx}{d \tau} &= 1 - x - \nu x f(x-y) - A\delta x \\
\frac{dy}{d \tau} &= \delta (R-y) - \nu y f(x-y) \,,
\end{aligned}
```

where ``x`` and ``y`` denote non-dimensional gradients of temperature and salinity, respectively, between the equatorial and the polar box.

The original Stommel 1961 model (see [`stommel`](@ref)), defined in terms of the parameters ``\\delta``, ``\\lambda``, and ``R``, is obtained by setting ``\\tau = ct``, ``\\nu = 1/\\lambda``, ``f(x-y) = |x-y|``, and ``A=0``.

Cessi's 1994 version of the model ([`cessi`](@ref)), written in terms of the parameters ``\\alpha``, ``\\mu^2``, and ``p``, follows from the above equations by setting ``\\tau = t/t_r``, ``\\delta = 2/\\alpha``, ``\\nu=\\mu^2``, ``f(x-y) = (x-y)^2/k``, ``R=pH/(2S_0)``, and ``A=1``.

```@docs
stommel(u, p, t; kwargs...)
cessi(u, p, t)
```

#### Interhemispheric 3-box model
```@docs
rooth_smooth(u, p, t; kwargs...)
```
