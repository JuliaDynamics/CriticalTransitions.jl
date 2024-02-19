"""
Script to generate an ensemble of noise-induced transitions in the fitzhugh_nagumo system
Author: Reyk Börner (reyk.boerner@reading.ac.uk)
Date: 3 Oct 2022

Note: To enable multi-threading, run the file with "julia --threads N sample_transitions.jl",
where N is to be replaced with the desired number of cores.

Note: To access the sampling data in Julia, use the JLD2 package and do 
data = jldopen("filepath/filename.jld2")
"""

# Load modules
include("../src/CriticalTransitions.jl")
using .CriticalTransitions, Printf


###########################################################
# SETTINGS ################################################
###########################################################

# fitzhugh_nagumo parameters
σ = 0.15                # noise intensity
ϵ = 1.0                 # time scale parameter

β, α, γ, κ, Ι = 3., 1., 1., 1., 0.
pf = [ϵ, β, α, γ, κ, Ι]

# Noise settings
Σ = [1 0; 0 1]          # covariance matrix
process = "WhiteGauss"  # noise process

# Transition settings
AtoB = true             # if false, transition B->A
rad_i = 0.03            # ball radius around initial point
rad_f = 0.1             # ball radius around final point

# Run settings
N = 500                	# number of samples
tmax = Int64(1e3)       # maximum simulation time
dt = 0.01               # time step
solver = EM()		# SDEProblem solver

# I/O settings
filetext = "fhn_transAB_eps$(@sprintf "%.2e" ϵ)_noise$(@sprintf "%.2e" σ)_eye_$(N)"
save_path = "../data/"

###########################################################
# End of settings. ########################################
###########################################################

# Instantiate system
sys = StochSystem(fitzhugh_nagumo, pf, 2, σ, idfunc, nothing, Σ, process)

# Get fixed points
pts, eigs, stab = fixedpoints(sys, [-10,-10], [10,10])
A, B = pts[stab]

# A->B or B->A?
if AtoB
    A, B = B, A
end

# I/O info
info0 = "$(N) transition samples of $(sys.f) system with noise intensity $(sys.σ)."
info1 = "Start/end conditions:\n Initial state: $(A), ball radius $(rad_i)\n Final state: $(B), ball radius $(rad_f)"
info2 = "Run details:\n time step dt=$(dt), maximum duration tmax=$(tmax), solver=$(string(solver))"
info3 = "Dataset info:\n paths: transition path arrays of dimension (state × time), where a state is a vector [u,v]\n times: arrays of time values of the corresponding paths"

println(info0*"\n... "*info1*"\n... "*info2)
println("... Saving data to yymmdd_$(save_path*filetext).jld2.")

# Create data file
file = make_jld2(filetext, save_path)
write(file, "system_info", sys_string(sys, verbose=true))
write(file, "run_info", info0*"\n"*info1*"\n"*info2)
write(file, "data_info", info3)

# Call transitions function
precompile(transitions, (StochSystem, State, State, Int64,))
samples, times, idx = transitions(sys, A, B, N;
    rad_i=rad_i, rad_f=rad_f, dt=dt, tmax=tmax, solver=solver,
    savefile=file, Nmax=1000)

# Save which sample numbers transitioned
write(file, "sample_idxs", string(sort(idx)))

# Close file
println("... Done! Closing file.")
println("... Summary: In $(length(idx)) of $(N) samples, a transition occurred.")
close(file)

# End of script ###########################################
