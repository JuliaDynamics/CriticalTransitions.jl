using DynamicalSystems

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

function intervals_to_box(bmin::Vector, bmax::Vector, dim::Int64)
    # Generates a box from specifying the interval limits
    intervals = []
    for i in 1:dim
        push!(intervals, bmin[i]..bmax[i])
    end
    box = intervals[1]
    for i in 2:dim
        box = box × intervals[i]
    end
    box
end;