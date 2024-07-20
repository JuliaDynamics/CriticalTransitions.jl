module CTLibrary

using CriticalTransitions

include("fitzhughnagumo.jl")
include("truscottbrindley_mod.jl")
include("truscottbrindley_orig.jl")
include("truscottbrindley_orig1.jl")
include("rooth.jl")
include("stommel.jl")
include("rivals.jl")

export fitzhugh_nagumo

end
