"""
CriticalTransitions function test file
Run this file to see whether any functions throw an error
"""

include("../src/CriticalTransitions.jl")
using .CriticalTransitions

# Set up an example system
pf = [1., 3., 1., 1., 1., 0.]
sys = StochSystem(FitzHughNagumo, pf, 2, 0.3, idfunc, nothing, I(2), "WhiteGauss")
sys! = StochSystem(FitzHughNagumo!, pf, 2, 0.3, idfunc!, nothing, I(2), "WhiteGauss")
# Check functions
eq = equilib(sys, [1.,1.])
fp = fixedpoints(sys, [-2,-2], [2,2])
sm = simulate(sys, [1.,1.])
rx = relax(sys, [1.,1.])
tr = transition(sys, fp[1][1], fp[1][3])
tt = transitions(sys, fp[1][1], fp[1][3], 2, rad_i=0.1, rad_f=0.1, tmax=1e3)
ds = tocds(sys)
st = sys_string(sys)
io = sys_info(sys)

println(EM())
println(I(2))