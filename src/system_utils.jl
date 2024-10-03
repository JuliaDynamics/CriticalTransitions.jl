"""
$(TYPEDSIGNATURES)

Returns the deterministic drift ``f(x)`` of the CoupledSDEs `sys` at state `x`. For time-
dependent systems, the time can be specified as a keyword argument `t` (by default `t=0`).
"""
function drift(sys::CoupledSDEs{IIP}, x; t=0) where {IIP}
    f = dynamic_rule(sys)
    if IIP
        dx = similar(x)
        f(dx, x, sys.p0, t)
        return dx
    else
        return f(x, sys.p0, t)
    end
end

"""
$(TYPEDSIGNATURES)

Computes the divergence of the drift field ``f(x)`` at state `x`. For time-
dependent systems, the time can be specified as a keyword argument `t` (by default `t=0`).
"""
function div_drift(sys::CoupledSDEs, x; t=0)
    b(x) = drift(sys, x; t)
    return tr(ForwardDiff.jacobian(b, x))
end;