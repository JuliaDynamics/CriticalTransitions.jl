"""
CriticalTransitions function test file
Run this file to see whether any functions throw an error
"""

include("../src/CriticalTransitions.jl")
using .CriticalTransitions

# Set up an example system
p = [1., 3., 1., 1., 1., 0.];
sys = StochSystem(fitzhugh_nagumo, p, zeros(2), 0.2, idfunc, nothing, I(2), "WhiteGauss")
sys! = StochSystem(fitzhugh_nagumo!, p, zeros(2), 0.3, idfunc!, nothing, I(2), "WhiteGauss")
# Check functions
ds = CoupledODEs(sys)
StochSystem(ds)
eq = equilib(sys, [1.,1.])
fp, eig, stab = fixedpoints(sys, [-2,-2], [2,2])
ba = basins(sys, [0.,0], [0.,1], [1.,0], intervals_to_box([-2,-2],[2,2]), bstep=[0.5,0.5], ϵ_mapper=0.001, Ttr=100)
bb = basinboundary(ba)
sm = simulate(sys, [1.,1.])
rx = relax(sys, [1.,1.])
tr = transition(sys, fp[stab][1], fp[stab][2])
tt = transitions(sys, fp[stab][1], fp[stab][2], 50, rad_i=0.1, rad_f=0.1, tmax=1e3)
mm = min_action_method(sys, fp[stab][1], fp[stab][2], 50, 3, maxiter=5)
gm = geometric_min_action_method(sys, fp[stab][1], fp[stab][2], N=50, maxiter=5)
et = edgetracking(sys, [-1,0.5], [1,0.5], [fp[stab][1], fp[stab][2]], maxit=5)
dr = drift(sys, [1,1])
st = sys_string(sys)
io = sys_info(sys)
R, L = fp[stab]
initial = reduce(hcat, range(R,L, length=50))
lv = langevinmcmc(sys, initial; tmax=0.1)

println(EM())
println(I(2))