using Formatting, Dates, JLD2

include("StochSystem.jl")

function make_jld2(text::String, relpath::String="")
    """
    Creates/opens a .jld2 file with date stamp and specified filename.
        Format: "ddmmyy_text.jdl2"
        Relative path specified by relpath (ending with "/").
    """
    time = Dates.now()
    str = relpath*Dates.format(time, "yymmdd")*"_"*text*".jld2"
    jldopen(str, "a+")
end;

function sys_info(sys::StochSystem)
    # Prints StochSystem info in clean, readable format
printfmt("StochSystem
 - f (deterministic function):  {:}
 - pf (parameters of f):        {:}
 - g (noise function):          {:}
 - pg (parameters of g):        {:}
 - σ (noise intensity):         {:.3e}
 - Σ (covariance matrix):       {:}
 - process (noise process):     {:}", sys.f, sys.pf, sys.g, sys.pg, sys.σ, sys.Σ, sys.process)
end;

function sys_string(sys::StochSystem; verbose=true)
    # Returns StochSystem info as a String
    if verbose
        "StochSystem\n - f (deterministic function): $(sys.f)\n - pf (parameters of f): $(sys.pf)\n - g (noise function): $(sys.g)\n - pg (parameters of g): $(sys.pg)\n - σ (noise intensity): $(sys.σ)\n - Σ (covariance matrix): $(sys.Σ)\n - process (noise process): $(sys.process)"
    else
        string(sys)
    end
end;