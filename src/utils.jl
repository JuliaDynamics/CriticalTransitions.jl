"""
Utility functions for CriticalTransitions.jl
Functions in this file are independent of types/functions defined in CriticalTransitions
"""

function idfunc(u, p, t)
    # Identity function for noise function g (out-of-place)
    ones(length(u))
end;

function idfunc!(du, u, p, t)
    # Identity function for noise function g (in-place)
    du = ones(length(u))
end;

function is_iip(f::Function)
    # Asserts if f is in-place (true) or out-of-place (false)
    occursin("!", String(Symbol(f)))
end;

function intervals_to_box(bmin::Vector, bmax::Vector)
    # Generates a box from specifying the interval limits
    intervals = []
    dim = length(bmin)
    for i in 1:dim
        push!(intervals, bmin[i]..bmax[i])
    end
    box = intervals[1]
    for i in 2:dim
        box = box × intervals[i]
    end
    box
end;

function additive_idx!(du, u, p, t, idx)
    du[indx] .= 1.
end;

function additive_idx(u,p,t,idx)
    du = zeros(length(u))
    du[idx] .= 1.
    SVector{length(u)}(du)
end;

function multiplicative_idx!(du, u, p, t, idx)
    du[indx] .= u[idx]
end;

function multiplicative_idx(u, p, t, idx)
    du = zeros(length(u))
    du[idx] = u[idx]
    SVector{length(u)}(du)
end


function multiplicative_first!(du, u, p, t)
    du[1] = u[1];
end;

function multiplicative_first(u, p, t)
    
    du = zeros(length(u));
    du[1] = u[1];

    SVector{length(u)}(du)
end;

function additive_first!(du, u, p, t)
    du[1] = 1;
end;

function additive_first(u, p, t)
    
    du = zeros(length(u));
    du[1] = 1;

    SVector{length(u)}(du)
end;