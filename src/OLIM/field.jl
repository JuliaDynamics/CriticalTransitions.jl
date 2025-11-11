
# // choose the vector field b
# const char chfield='l';
# // Vector fields with asymptotically stable equilibrium
# // chfield = 'l': b = [-2, -alin; 2*alin, -1][ x; y]; (in the Matlab notations),
# // analytic solution  U = 2x^2 + y^2
# const double alin = 10.0;
# // chfield = 'q' ->  Maier-Stein model;
# // chfield = 'r' -> FitzHugh-Nagumo model with a = 1.2

# // Vector fields with stable limit cycle
# // chfield = 'b' ->  Brusselator;
# // CHANGE const char *ficurve_name = "circle.txt"; to
# // const char *ficurve_name = "icurve_brusselator.txt";
# // chfield = 'c' -> analytic solution U = (1/2)(r^2-1)^2, l = (y,-x)

const alin = 10.0;
function myfield(x::Vector{Float64}, chfield::Char)::Vector{Float64}
    v = zeros(2)

    if chfield == 'b'
        v[1] = 1.0 + x[1]^2 * x[2] - 4.0 * x[1]
        v[2] = 3 * x[1] - x[1]^2 * x[2]
    elseif chfield == 'c'
        aux = 1.0 - x[1]^2 - x[2]^2
        v[1] = x[2] + x[1] * aux
        v[2] = -x[1] + x[2] * aux
    elseif chfield == 'q'
        v[1] = x[1] - x[1]^3 - 10.0 * x[1] * x[2]^2
        v[2] = -(1.0 + x[1]^2) * x[2]
    elseif chfield == 'r'
        v[1] = x[1] - x[1]^3/3.0 - x[2]
        v[2] = x[1] + 1.2
    elseif chfield == 'l'
        v[1] = -2 * x[1] - alin * x[2]
        v[2] = 2 * alin * x[1] - x[2]
    else
        error("Invalid chfield: $chfield")
    end

    return v
end

function exact_solution(x::Vector{Float64}, chfield::Char)::Float64
    if chfield == 'l'
        return 2.0 * x[1]^2 + x[2]^2
    elseif chfield == 'c'
        return 0.5 * (x[1]^2 + x[2]^2 - 1.0)^2
    else
        return 0.0
    end
end

function param(chfield::Char)
    # Define parameters based on the field type
    x_ipoint = Dict(
        'q' => [-1.0, 0.0],
        'r' => [-1.2, -1.2*(1.0-1.2*1.2/3.0)],
        'l' => [0.0, 0.0]
    )

    bounds = Dict(
        'q' => (-2.0, 2.0, -2.0, 2.0),
        'r' => (-2.5, 2.5, -2.5, 2.5),
        'l' => (-1.0, 1.0, -1.0, 1.0),
        'c' => (-2.0, 2.0, -2.0, 2.0),
        'b' => (0.0, 7.0, 0.0, 7.0)
    )

    if !haskey(bounds, chfield)
        error("Invalid chfield: $chfield")
    end

    XMIN, XMAX, YMIN, YMAX = bounds[chfield]

    println("in param()")

    # Grid parameters
    hx = (XMAX - XMIN)/(NX-1)
    hy = (YMAX - YMIN)/(NY-1)
    h = sqrt(hx^2 + hy^2)

    # Initialize arrays
    ms = zeros(Int, NX * NY)
    g = fill(Inf, NX * NY)
    solinfo = fill('n', NX * NY)

    # Compute exact solution if available
    Uexact = zeros(NX * NY)
    if chfield in ['l', 'c']
        for j in 0:NY-1
            for i in 0:NX-1
                ind = i + 1 + NX * j
                x = [XMIN + i*hx, YMIN + j*hy]
                Uexact[ind] = exact_solution(x, chfield)
            end
        end
    end

    return hx, hy, h, ms, g, solinfo, Uexact, x_ipoint
end
