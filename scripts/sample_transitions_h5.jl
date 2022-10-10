"""
Script to generate an ensemble of noise-induced transitions in the FitzHughNagumo system
Author: Reyk Börner (reyk.boerner@reading.ac.uk)
Date: 3 Oct 2022

Note: To enable multi-threading, run the file with "julia --threads N sample_transitions.jl",
where N is to be replaced with the desired number of cores.

Note: This script saves the data in HDF5 format.
To access the sampling data in Julia, use the HDF5 package
and do data = h5open("filepath/filename.h5", "r").
"""

# Load modules
using CriticalTransitions, Printf, HDF5, Dates


###########################################################
# SETTINGS ################################################
###########################################################

# FitzHughNagumo parameters
σ = 0.22               # noise intensity
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
N = 25              	# number of transition samples
Nmax = 100              # maximum number of attempts
tmax = Int64(1e4)       # maximum simulation time per attempt
dt = 0.02               # time step
solver = EM()		# SDEProblem solver

# I/O settings
filetext = "fhn_transAB_eps$(@sprintf "%.2e" ϵ)_noise$(@sprintf "%.2e" σ)_eye_$(N)"
save_path = "../data/"

###########################################################
# End of settings. ########################################
###########################################################

# Instantiate system
sys = StochSystem(FitzHughNagumo, pf, 2, σ, idfunc, nothing, Σ, process)

# Fixed points
A, B = [-sqrt(2/3), -sqrt(2/27)], [sqrt(2/3), sqrt(2/27)]

# A->B or B->A?
if AtoB
    A, B = B, A
end

# I/O info
info0 = "$(N) transition samples of $(sys.f) system with noise intensity $(sys.σ)."
info1 = "Start/end conditions:\n Initial state: $(A), ball radius $(rad_i)\n Final state: $(B), ball radius $(rad_f)"
info2 = "Run details:\n time step dt=$(dt), maximum duration tmax=$(tmax), solver=$(string(solver))"
info3 = "Dataset info:\n paths: transition path arrays of dimension (state × time), where a state is a vector [u,v]\n times: arrays of time values of the corresponding paths"

println(info0*"\n--> "*info1*"\n--> "*info2)
println("--> Saving data to $(save_path).")

# Create data file
#file = make_h5(filetext, save_path)
file = make_h5(filetext)

create_group(file, "paths")
create_group(file, "times")
paths = file["paths"]
times = file["times"]
attributes(paths)["data dimensions"] = "[coordinate (rows) × time (columns)], coordinates = [u,v], arbitrary units"
attributes(times)["data dimensions"] = "[time (since beginning of simulation)], arbitrary units"
attributes(file)["system_info"] = sys_string(sys, verbose=true)
attributes(file)["data_info"] = info3

# Call transitions function
tstart = now()
samples, times, idx = transitions(sys, A, B, N;
    rad_i=rad_i, rad_f=rad_f, dt=dt, tmax=tmax, solver=solver,
    savefile=file, Nmax=Nmax)

# Finalize file
attributes(file)["run_info"] = info0*"\n"*info1*"\n"*info2*"\nRuntime: $(now()-tstart)"
write(file, "path_numbers", sort(idx))

# Close file
println("... Done! Closing file. Run took $(now()-tstart).")
close(file)
println("... Summary: $(length(idx)) samples transitioned.")


# End of script ###########################################