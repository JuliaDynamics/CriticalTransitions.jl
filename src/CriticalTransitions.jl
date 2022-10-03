module CriticalTransitions

include("StochSystem.jl")
include("utils.jl")
include("stability.jl")
include("simulation.jl")

include("../systems/fitzhughnagumo.jl")

export StochSystem
export equilib, fixedpoints
export simulate, relax, transition, transitions
export tocds, make_jld2

export FitzHughNagumo, FitzHughNagumo!
export idfunc, idfunc!

end # module CriticalTransitions