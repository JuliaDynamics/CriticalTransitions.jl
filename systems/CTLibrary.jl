module CTLibrary

using CriticalTransitions: CriticalTransitions
using IntervalArithmetic: IntervalArithmetic, interval

include("fitzhughnagumo.jl")
include("truscottbrindley_mod.jl")
include("truscottbrindley_orig.jl")
include("truscottbrindley_orig1.jl")
include("rooth.jl")
include("stommel.jl")
include("rivals.jl")
include("ou_multiplicative.jl")

export fitzhugh_nagumo!, fitzhugh_nagumo, stommel, rivals!, rivals, cessi, rooth_smooth
export ou_multiplicative_1d, linear_offdiag_2d_sde
#modifiedtruscottbrindleywithdimensions!, modifiedtruscottbrindleywithdimensions
#originaltruscottbrindley!, originaltruscottbrindley
#rampedoriginaltruscottbrindley!, rampedoriginaltruscottbrindley

end
