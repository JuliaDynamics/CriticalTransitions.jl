# test_include.jl - stub environment to parse RateSystem.jl
module DynamicalSystemsBase
export current_parameters, current_parameter, dynamic_rule, initial_time, current_state, set_parameter!, set_parameters!, initial_state, initial_parameters, current_time, successful_step, set_state!, trajectory

function initial_time(ds)
	return 0.0
end

function dynamic_rule(ds)
	return (du, u, p, t) -> du
end

function current_parameter(ds, k)
	return 0.0
end

function current_parameters(ds)
	return Dict{Any,Any}()
end

function current_state(ds)
	return zeros(0)
end

function set_parameter!(ds, k, v)
	return ds
end

function set_parameters!(ds, params)
	return ds
end

function initial_state(ds)
	return zeros(0)
end

function initial_parameters(ds)
	return Dict{Any,Any}()
end

function current_time(ds)
	return 0.0
end

function successful_step(ds)
	return true
end

function set_state!(ds, u)
	return ds
end

function isdeterministic(ds)
	return true
end

function trajectory(ds, args...; kw...)
	return (Array{Float64,2}(undef,0,0), Float64[])
end

end
module SciMLBase
export step!, isinplace

function step!(ds, args...)
	return ds
end

function isinplace(ds)
	return false
end

end
module StateSpaceSets
export dimension

function dimension(ds)
	return 0
end

end

abstract type ContinuousTimeDynamicalSystem end
struct CoupledODEs end
CoupledODEs(args...; kwargs...) = CoupledODEs()

function initial_time(ds)
	return 0.0
end

function dynamic_rule(ds)
	return (du, u, p, t) -> du
end

function current_parameter(ds, k)
	return 0.0
end

function current_parameters(ds)
	return Dict{Any,Any}()
end

function current_state(ds)
	return zeros(0)
end

# Provide placeholders for interpolated docstring variables used in sources
const TYPEDSIGNATURES = ""

include(joinpath(@__DIR__, "src", "r_tipping", "RateSystem.jl"))
println("Included RateSystem.jl")
