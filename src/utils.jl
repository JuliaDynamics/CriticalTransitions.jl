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