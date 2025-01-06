function basins end
function basboundary end
function basinboundary end
function intervals_to_box end

function get_boundary end
function find_boundary end
function huniform end
function dellipse end
function dunion end
function ddiff end
function distmesh2D end
function get_ellipse end
function reparametrization end

# ## Method error handling
# We also inject a method error handler, which
# prints a suggestion if the Proj extension is not loaded.
function _basin_error_hinter(func)
    return function _dump_error_hinter(io, exc, argtypes, kwargs)
        if isnothing(Base.get_extension(CriticalTransitions, :CoupledSDEsBasin)) &&
            exc.f == func
            print(
                io,
                "\n\nThe `$(string(func))` method requires the ChaosTools.jl and Attractors.jl package to be explicitly loaded.\n",
            )
            print(io, "You can do this by simply typing ")
            printstyled(io, "using ChaosTools, Attractors"; color=:cyan, bold=true)
            println(
                io,
                " in your REPL, \nor otherwise loading ChaosTools.jl and Attractors.jl via using or import.",
            )
        else # this is a more general error
            nothing
        end
    end
end
