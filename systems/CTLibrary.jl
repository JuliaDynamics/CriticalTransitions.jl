module CTLibrary

using CriticalTransitions

include("fitzhughnagumo.jl")
include("truscottbrindley_mod.jl")
include("truscottbrindley_orig.jl")
include("truscottbrindley_orig1.jl")
include("rooth.jl")
include("stommel.jl")
include("rivals.jl")

export fitzhughnagumo, fitzhughnagumo!, stommel, rivals!, rival, cessi, rooth_smooth
modifiedtruscottbrindleywithdimensions!, modifiedtruscottbrindleywithdimensions
originaltruscottbrindley!, originaltruscottbrindley
rampedoriginaltruscottbrindley!, rampedoriginaltruscottbrindley

end
