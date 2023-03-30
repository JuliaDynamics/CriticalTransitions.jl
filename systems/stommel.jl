#include("../../src/utils.jl")

"""
    stommel(u, p, t; kwargs...)
Stommel's hemispheric 2-box model of the Thermohaline Circulation
([Stommel 1961](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.2153-3490.1961.tb00079.x)).

## Parameters 
`p = [[delta, mu, R]]`

## Keyword arguments
* `smooth = 1e6`: Approximates ``|x| \\approx x \\tanh(sx)``, where ``s =`` `smooth`.
  If `smooth=0`, the exact absolute value is taken, making the equations non-smooth.
* `flow_law`: Relation between density difference and flow strength.
  Options:
  * `"abs"`: Original Stommel model [Stommel 1961]
  * `"diffu_abs"`: Absolute flow law Q2 in eq. (2.4) of Cessi 1994[^1]
  * `"diffu_sqr"`: Absolute flow law Q3 in eq. (2.4) of Cessi 1994[^1]

See also [`cessi`](@Ref).

[^1]: [Cessi (1994)](https://journals.ametsoc.org/view/journals/phoc/24/9/1520-0485_1994_024_1911_asbmos_2_0_co_2.xml)
"""
function stommel(u, p, t; smooth=1e6, flow_law="abs")
    x, y = u
    delta, mu, R = p[1]

    if flow_law == "abs"
        # Original Stommel model
        q = (smooth>0) ? smoothabs(x-y, smooth) : abs(x-y)
        diffu = 0
    elseif flow_law == "diffu_abs"
        # Absolute flow law Q2 in eq. (2.4) of Cessi 1994
        q = (smooth>0) ? smoothabs(x-y, smooth) : abs(x-y)
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