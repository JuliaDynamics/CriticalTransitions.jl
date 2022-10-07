include("../StochSystem.jl")
include("../utils.jl")

function gauss(sys::StochSystem)
    # Returns a Wiener process for given covariance matrix and dimension of a StochSystem
    if is_iip(sys.f)
        W = CorrelatedWienerProcess!(sys.Σ, 0.0, zeros(sys.dim))
    else
        W = CorrelatedWienerProcess(sys.Σ, 0.0, zeros(sys.dim))
    end
    W
end;