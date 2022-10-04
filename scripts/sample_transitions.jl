"""
Script to generate an ensemble of noise-induced transitions in the FitzHughNagumo system
Author: Reyk Börner
Date: 3 Oct 2022

Note: To enable multi-threading, run the file with "julia --nthreads N sample_transitions.jl",
where N is to be replaced with the number of cores to be used.
"""
# Load modules
using DifferentialEquations, Printf
include("../src/CriticalTransitions.jl")
using .CriticalTransitions

###########################################################
# SETTINGS ################################################
###########################################################

# FitzHughNagumo parameters
σ = 0.18                # noise intensity
ϵ = 1.0                 # time scale parameter

β, α, γ, κ, Ι = 3., 1., 1., 1., 0.
pf = [ϵ, β, α, γ, κ, Ι]

# Noise settings
Σ = [1 0; 0 1]          # covariance matrix
process = "WhiteGauss"  # noise process

# Transition settings
AtoB = true             # if false, transition B->A
rad_i = 0.05            # ball radius around initial point
rad_f = 0.1             # ball radius around final point

# Run settings
N = 100                 # number of samples
tmax = Int64(1e3)       # maximum simulation time
dt = 0.01               # time step
solver = EM()			# SDEProblem solver

# I/O settings
filetext = "fhn_transAB_eps$(@sprintf "%.2e" ϵ)_noise$(@sprintf "%.2e" σ)_eye_$(N)"
save_path = "../data/"

###########################################################
# End of settings. ########################################
###########################################################

# Instantiate system
sys = StochSystem(FitzHughNagumo, pf, 2, σ, idfunc, nothing, Σ, process)

# Get fixed points
pts, eigs, stab = fixedpoints(sys, [-10,-10], [10,10])
A, B = pts[stab]

# A->B or B->A?
if AtoB
    A, B = B, A
end

# Store info
info0 = "$(N) transition samples of $(sys.f) system with noise intensity $(sys.σ)."
info1 = "Start/end conditions:\n Initial state: $(A), ball radius $(rad_i)\n Final state: $(B), ball radius $(rad_f)"
info2 = "Run details:\n time step dt=$(dt), maximum duration tmax=$(tmax), solver=$(string(solver))"

println(info0*"\n... "*info1*"\n... "*info2)
println("... saving data to $(save_path*filetext).jld2.")

# Create data file
file = make_jld2(filetext, save_path)
write(file, "system_info", sys_string(sys, verbose=true))
write(file, "run_info", info0*"\n"*info1*"\n"*info2)

# Call transitions function
samples, times, idx = transitions(sys, A, B, N;
    rad_i=rad_i, rad_f=rad_f, dt=dt, tmax=tmax, solver=solver,
    savefile=file)

# Save which sample numbers transitioned
write(file, "sample_idxs", sort(idx))

# Close file
println("... Done! Closing file...")
close(file)

# End of script ###########################################