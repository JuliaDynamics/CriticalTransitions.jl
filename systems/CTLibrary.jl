module CTLibrary

using CriticalTransitions: CriticalTransitions, smoothabs
using IntervalArithmetic: interval
using StaticArrays: SA, SVector

include("fitzhughnagumo.jl")
include("truscottbrindley_mod.jl")
include("truscottbrindley_orig.jl")
include("truscottbrindley_orig1.jl")
include("rooth.jl")
include("stommel.jl")
include("rivals.jl")

export fitzhugh_nagumo!, fitzhugh_nagumo, stommel, rivals!, rivals, cessi, rooth_smooth
modifiedtruscottbrindleywithdimensions!, modifiedtruscottbrindleywithdimensions
originaltruscottbrindley!, originaltruscottbrindley
rampedoriginaltruscottbrindley!, rampedoriginaltruscottbrindley

end
