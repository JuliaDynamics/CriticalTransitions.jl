"""
$(TYPEDSIGNATURES)

Creates/opens a `.jld2` file with filename of the format "ddmmyy_`text`.jld2". Relative file path specified by `relpath` (must end with `/`).

See [`make_h5`](@ref) for generating a `.h5` file.
"""
function make_jld2(text::String, relpath::String="")
    time = Dates.now()
    str = relpath*Dates.format(time, "yymmdd")*"_"*text*".jld2"
    jldopen(str, "a+")
end;

"""
$(TYPEDSIGNATURES)

Creates/opens a `.h5` file with filename of the format "ddmmyy_`text`.h5". Relative file path specified by `relpath` (must end with `/`).

See [`make_jld2`](@ref) for generating a `.jld2` file.
"""
function make_h5(text::String, relpath::String="")
    time = Dates.now()
    str = relpath*Dates.format(time, "yymmdd")*"_"*text*".h5"
    h5open(str, "cw")
end;
