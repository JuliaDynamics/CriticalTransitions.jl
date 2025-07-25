# Setting up a `RateSystem`
Let us explore an example of how to construct a `RateSystem` from an autonomous dynamical system (e.g. a `CoupledODEs`) and a time-dependent forcing protocol called `RateProtocol`.

```@example RateSystem
using CriticalTransitions
using CairoMakie
using CairoMakie.Makie.MathTeXEngine: get_font
font = (;
    regular=get_font(:regular), bold=get_font(:bold),
    italic=get_font(:italic), bold_italic=get_font(:bolditalic)
);
```

## Prototypical model for R-tipping, with critical rate r = 4/3

We first consider the following simple one-dimensional autonomous system with one attractor, given by the ordinary differential equation:
```math
\begin{aligned}
    \dot{x} &= (x+\lambda)^2 - 1
\end{aligned}
```
The parameter ``\lambda`` shifts the location of the extrema of the drift field. 
We implement this system as follows:

```@example RateSystem
function f(u,p,t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end

lambda = 0.0 
p = [lambda]
x0 = [-1.]
auto_sys = CoupledODEs(f,x0,p)
```

## Non-autonomous case

Now, we want to explore a non-autonomous version of the system. 
We consider a setting where in the past and in the future the system is autnonomous and in between there is a non-autonomous period ``[t_start, t_end]`` with a time-dependent parameter ramping given by the function ``\lambda(rt)``. Choosing different values of the parameter ``r`` allows us to vary the speed of the parameter ramping.

We start by defining the function ``\lambda(t)``:
```@example RateSystem
function λ(p,t)
    λ_max = p[1]
    lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
    return SVector{1}(lambda)
end

λ_max = 3.
p_lambda = [λ_max] # parameter of the function lambda
```

Now, we define the RateProtocol that describes the non-autonomous period:

```@example RateSystem
r = 4/3-0.02   # r just below critical rate
t_start = -Inf # start time of non-autonomous part
t_end = Inf    # end time of non-autonomous part

rp = CriticalTransitions.RateProtocol(λ,p_lambda,r,t_start,t_end)
```

Now, we set up the combined system with autonomous past and future and non-autonomous ramping in between:

```@example RateSystem
t0 = -10.      # initial time of the system
nonauto_sys = RateSystem(auto_sys,rp,t0)
```

We can compute trajectories of this new system in the familiar way:
```@example RateSystem
T = 20.        # final simulation time
auto_traj = trajectory(auto_sys,T,x0)
nonauto_traj = trajectory(nonauto_sys,T,x0)
```

We plot the two trajectories

```@example RateSystem
fig = Figure(); axs = Axis(fig[1,1])
lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}")
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}")
axislegend(axs,position=:rc,labelsize=10)
fig
```

-----
Author: Raphael Roemer

Date: 30 Jun 2025
