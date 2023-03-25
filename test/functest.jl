"""
CriticalTransitions function test file
Run this file to see whether any functions throw an error
"""

include("../src/CriticalTransitions.jl")
using .CriticalTransitions

# Set up an example system
p = [1., 3., 1., 1., 1., 0.];
sys = StochSystem(FitzHughNagumo, p, 2, 0.2, idfunc, nothing, I(2), "WhiteGauss")
sys! = StochSystem(FitzHughNagumo!, p, 2, 0.3, idfunc!, nothing, I(2), "WhiteGauss")
# Check functions
eq = equilib(sys, [1.,1.])
fp, eig, stab = fixedpoints(sys, [-2,-2], [2,2])
ba = basins(sys, [0.,0], [0.,1], [1.,0], intervals_to_box([-2,-2],[2,2]), bstep=[0.5,0.5])
bb = basinboundary(ba)
sm = simulate(sys, [1.,1.])
rx = relax(sys, [1.,1.])
tr = transition(sys, fp[stab][1], fp[stab][2])
tt = transitions(sys, fp[stab][1], fp[stab][2], 50, rad_i=0.1, rad_f=0.1, tmax=1e3)
mm = mam(sys, fp[stab][1], fp[stab][2], 50, 3, maxiter=5)
gm = gmam(sys, fp[stab][1], fp[stab][2], N=50, maxiter=5)
et = edgetracking(sys, [-1,0.5], [1,0.5], [fp[stab][1], fp[stab][2]], maxit=5)
dr = drift(sys, [1,1])
ds = to_cds(sys)
st = sys_string(sys)
io = sys_info(sys)

println(EM())
println(I(2))