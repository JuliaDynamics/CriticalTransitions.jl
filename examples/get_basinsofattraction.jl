"""
Script to compute and save the basin boundary of the fitzhugh_nagumo system
for different values of time scale separation
"""

# Load packages
include("../src/CriticalTransitions.jl")
using .CriticalTransitions, HDF5

# Define systems
eps = [0.01, 0.1, 1., 10.];
systems = Dict();
for i in eps
    systems[i] = StochSystem(fitzhugh_nagumo, [i, 3, 1, 1, 1, 0], 2, 0.1)
end;

# Set up domain
A, B, C = [0.,0.], [1.,0.], [0.,1.];
box = intervals_to_box([-2,-1.5], [2,1.5]);
step = 0.02;

# Compute basins of attraction
function iterate_boas(p)
    _boas = Dict();
    Threads.@threads for i in p 
        _boas[i] = basins(systems[i], A, B, C, box; bstep=[step, step])
    end
    _boas
end
boas = iterate_boas(eps)

# Compute basin boundary
bnds = Dict();
for i in eps
    bnds[i] = basinboundary(boas[i])
end;

file = h5open("fhn_basinboundaries_res$(step).h5", "cw");
attributes(file)["info"] = "Basins boundary of fitzhugh_nagumo system for different time scale parameter values epsilon";
attributes(file)["data_dimensions"] = "(dim × N), where dim = system dimension and N = number of points of curve. System coordinates = [u (row 1), v (row 2)]"
for i in eps
    write(file, "boundary_eps=$(i)", bnds[i])
    attributes(file)["systeminfo_eps=$(i)"] = sys_string(systems[i])
end;

close(file);
