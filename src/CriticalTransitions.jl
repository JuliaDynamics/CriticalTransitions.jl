module CriticalTransitions

include("StochSystem.jl")
include("utils.jl")
include("stability.jl")
include("simulation.jl")

include("../systems/fitzhughnagumo.jl")

export StochSystem
export equilib, fixedpoints
export simulate, relax, transition
export tocds

export FitzHughNagumo, FitzHughNagumo!

end # module CriticalTransitions