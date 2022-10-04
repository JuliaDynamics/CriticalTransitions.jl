module CriticalTransitions

include("StochSystem.jl")
include("utils.jl")
include("io.jl")
include("noise.jl")
include("stability.jl")
include("simulation.jl")

include("../systems/fitzhughnagumo.jl")

export StochSystem
export equilib, fixedpoints
export simulate, relax, transition, transitions
export tocds, make_jld2, sys_string, sys_info

export FitzHughNagumo, FitzHughNagumo!
export idfunc, idfunc!

end # module CriticalTransitions