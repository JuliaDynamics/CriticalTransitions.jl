#include("../StochSystem.jl")

"""
    make_jld2(text::String, relpath::String="")
Creates/opens a `.jld2` file with filename of the format "ddmmyy_`text`.jld2". Relative file path specified by `relpath` (must end with `/`).

See [`make_h5`](@ref) for generating a `.h5` file.
"""
function make_jld2(text::String, relpath::String="")
    time = Dates.now()
    str = relpath*Dates.format(time, "yymmdd")*"_"*text*".jld2"
    jldopen(str, "a+")
end;

"""
    make_h5(text::String, relpath::String="")
Creates/opens a `.h5` file with filename of the format "ddmmyy_`text`.h5". Relative file path specified by `relpath` (must end with `/`).

See [`make_jld2`](@ref) for generating a `.jld2` file.
"""
function make_h5(text::String, relpath::String="")
    time = Dates.now()
    str = relpath*Dates.format(time, "yymmdd")*"_"*text*".h5"
    h5open(str, "cw")
end;

"""
    sys_info(sys::StochSystem)
Prints StochSystem info of `sys` in structured format.
"""
function sys_info(sys::StochSystem)
printfmt("{:}-dimensional stochastic dynamical system
 - f (deterministic function):  {:}
 - pf (parameters of f):        {:}
 - g (noise function):          {:}
 - pg (parameters of g):        {:}
 - σ (noise intensity):         {:.3e}
 - Σ (covariance matrix):       {:}
 - process (noise process):     {:}", length(sys.u), sys.f, sys.pf, sys.g, sys.pg, sys.σ, sys.Σ,
 sys.process)
end;

"""
    sys_string(sys::StochSystem; verbose=true)
Returns StochSystem info of `sys` as a string.

## Keyword arguments
`verbose`: if true, the string includes additional descriptions. 
"""
function sys_string(sys::StochSystem; verbose=true)
    if verbose
        "StochSystem\n - f (deterministic function): $(sys.f)\n - pf (parameters of f): $(sys.pf)\n - g (noise function): $(sys.g)\n - pg (parameters of g): $(sys.pg)\n - σ (noise intensity): $(sys.σ)\n - Σ (covariance matrix): $(sys.Σ)\n - process (noise process): $(sys.process)"
    else
        string(sys)
    end
end;

# Pretty printing (based on https://github.com/JuliaDynamics/DynamicalSystemsBase.jl/blob/main/src/core/pretty_printing.jl)
Base.summary(sys::StochSystem) = "$(length(sys.u))-dimensional stochastic dynamical system"
function Base.show(io::IO, sys::StochSystem)
    print(io, summary(sys))
end

additional_details(::StochSystem) = []

function Base.show(io::IO, ::MIME"text/plain", sys::StochSystem)
    descriptors = [
        "f (deterministic function)" => sys.f,
        "pf (parameters of f)" => sys.pf,
        "g (noise function)" => sys.g,
        "pg (parameters of g)" => sys.pg,
        "σ (noise intensity)" => sys.σ,
        "Σ (covariance matrix)" => sys.Σ,
        "process (noise process)" => sys.process,
    ]
    padlen = maximum(length(d[1]) for d in descriptors) + 3

    println(io, summary(sys))
    for (desc, val) in descriptors
        println(io, rpad(" $(desc): ", padlen), val)
    end
end