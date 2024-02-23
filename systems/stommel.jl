"""
    stommel(u, p, t; kwargs...)
Stommel's hemispheric 2-box model of the Thermohaline Circulation[^Stommel1961].

The state vector `u = [x, y]` describes the non-dimensional temperature (`x`) and salinity (`y`) gradients
between the equatiorial and the polar box (details in the
[online docs](https://reykboerner.github.io/CriticalTransitions.jl/stable/man/systems/#Thermohaline-Circulation-box-models)).

The parameter vector `p = [[delta, mu, R]]` comprises the ratio ``\\delta`` of the saline and thermal
relaxation coefficients, the advective coefficient ``\\mu = 1/\\lambda``, and the constant ``R``
as described in [^Stommel 1961]. The original values in the bistable regime are `p = [[1/6, 5, 2]]`.

## Keyword arguments
* `smooth = 1e-6`: Approximates ``|x| \\approx x \\tanh(x/s)``, where ``s =`` `smooth`.
  If `smooth=0`, the exact absolute value is taken, making the equations non-smooth.
* `flow_law`: Relation between density difference and flow strength.
  Options:
  * `"abs"`: Original Stommel model [Stommel 1961]
  * `"diffu_abs"`: Absolute flow law Q2 in eq. (2.4) of [^Cessi1994]
  * `"diffu_sqr"`: Absolute flow law Q3 in eq. (2.4) of [^Cessi1994]

See also [`cessi`](@ref).

[^Stommel1961]:
    [Stommel (1961). Thermohaline Convection with two stable regimes of flow]
    (https://onlinelibrary.wiley.com/doi/abs/10.1111/j.2153-3490.1961.tb00079.x)
[^Cessi1994]:
    [Cessi (1994). A simple box model of stochastically forced thermohaline flow]
    (https://journals.ametsoc.org/view/journals/phoc/24/9/1520-0485_1994_024_1911_asbmos_2_0_co_2.xml)
"""
function stommel(u, p, t; smooth=1e-6, flow_law="abs")
    x, y = u
    delta, mu, R = p[1]

    if flow_law == "abs"
        # Original Stommel model
        q = (smooth>0) ? smoothabs(x-y, 1/smooth) : abs(x-y)
        diffu = 0
    elseif flow_law == "diffu_abs"
        # Absolute flow law Q2 in eq. (2.4) of Cessi 1994
        q = (smooth>0) ? smoothabs(x-y, 1/smooth) : abs(x-y)
        diffu = 1
    elseif flow_law == "diffu_sqr"
        # Quadratic flow law Q3 in eq. (2.4) of Cessi 1994
        q = (x-y)^2
        diffu = 1
    else
        @error("Invalid value of 'flow_law' kwarg. Options: 'abs', 'diffu_abs', 'diffu_sqr'.")
    end
    
    dx = 1 - x - mu*x*q - diffu*delta*x
    dy = delta*(R - y) - mu*y*q

    SA[dx, dy]
end

"""
    cessi(u, p, t)
Cessi's hemispheric 2-box model of the Thermohaline Circulation
([Cessi (1994)](https://journals.ametsoc.org/view/journals/phoc/24/9/1520-0485_1994_024_1911_asbmos_2_0_co_2.xml)).

## Parameters 
`p = [[alpha, mu2, pflux]]`

See also [`stommel`](@ref).
"""
function cessi(u, p, t)
    # Non-dimensional model parameters [Cessi 1994]
    alpha, mu2, pflux = p[1]
    
    # Convert to Stommel parameters
    delta = 2/alpha
    mu = 2/alpha*mu2
    R = pflux/2

    stommel(u, [[delta, mu, R]], t; smooth=0, flow_law="diffu_sqr")
end