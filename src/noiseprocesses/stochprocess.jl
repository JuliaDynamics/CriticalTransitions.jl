include("../StochSystem.jl")
include("../utils.jl")
include("gaussian.jl")

"""
    stochprocess(sys::StochSystem)

Translates the stochastic process specified in `sys` into the language required by the
`SDEProblem` of `DynamicalSystems.jl`.
"""
function stochprocess(sys::StochSystem)
    if sys.process == "WhiteGauss"
        if sys.Σ == I(sys.dim)
            return nothing
        else
            return gauss(sys)
        end
    else
        ArgumentError("Noise process not yet implemented.")
    end
end