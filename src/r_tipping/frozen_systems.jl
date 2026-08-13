# ! TODO write tests

"""
    frozen_system(rs::RateSystem, t)

Returns the autonomous version of the non-autonomous RateSystem `rs`
with parameters fixed at system time `t`.
"""
function frozen_system(rs::RateSystem, t)
    sys = rs.specs.unforced_system
    pnew = parameters(rs, t)
    try
        set_parameters!(sys, pnew)
    catch err
        for f in fieldnames(typeof(pnew))
            set_parameter!(sys, f, getfield(pnew, f))
        end
    end
    return sys
end

"""
    past_limit_system(rs::RateSystem, t)

Returns the autonomous dynamical system corresponding to the RateSystem `rs`
for times `t < forcing_start_time`. Equivalent to `rs.specs.unforced_system`.
"""
past_limit_system(rs::RateSystem) = rs.specs.unforced_system


"""
    future_limit_system(rs::RateSystem, t)

Returns the autonomous dynamical system corresponding to the RateSystem `rs`
for times `t > forcing_start_time + 2*forcing_duration`, i.e. after the time-dependent
forcing interval.
"""
future_limit_system(rs::RateSystem) = frozen_system(rs, Inf)