#include("../StochSystem.jl")
#include("../utils.jl")

"""
    gauss(sys::StochSystem)
Returns a Wiener process with dimension `sys.dim` and covariance matrix `sys.Σ`.

This function is based on the [`CorrelatedWienerProcess`](https://noise.sciml.ai/stable/noise_processes/#DiffEqNoiseProcess.CorrelatedWienerProcess) of [`DiffEqNoiseProcess.jl`](https://noise.sciml.ai/stable/), a component of `DifferentialEquations.jl`. The initial condition of the process is set to the zero vector at `t=0`.
"""
function gauss(sys::StochSystem)
    # Returns a Wiener process for given covariance matrix and dimension of a StochSystem
    if is_iip(sys.f)
        W = CorrelatedWienerProcess!(sys.Σ, 0.0, zeros(sys.dim))
    else
        W = CorrelatedWienerProcess(sys.Σ, 0.0, zeros(sys.dim))
    end
    W
end;

"""
    gauss(sys::RateSystem)
Returns a Wiener process with dimension `sys.dim` and covariance matrix `sys.Σ`.

This function is based on the [`CorrelatedWienerProcess`](https://noise.sciml.ai/stable/noise_processes/#DiffEqNoiseProcess.CorrelatedWienerProcess) of [`DiffEqNoiseProcess.jl`](https://noise.sciml.ai/stable/), a component of `DifferentialEquations.jl`. The initial condition of the process is set to the zero vector at `t=0`.
"""
function gauss(sys::RateSystem)
    # Returns a Wiener process for given covariance matrix and dimension of a StochSystem
    if is_iip(sys.f)
        W = CorrelatedWienerProcess!(sys.Σ, 0.0, zeros(sys.dim))
    else
        W = CorrelatedWienerProcess(sys.Σ, 0.0, zeros(sys.dim))
    end
    W
end;