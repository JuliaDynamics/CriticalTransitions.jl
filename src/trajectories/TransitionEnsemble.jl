"""
$(TYPEDEF)

Statistics of the ensemble of transition paths between two points in a state space.

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

"""
struct TransitionStatistics{T}
    "success rate of the transition process"
    success_rate::T
    "mean residence time of the transition process"
    residence_time::T
    "mean transition time of the transition process"
    transition_time::T
    "rareness of the transition process"
    rareness::T

    function TransitionStatistics(sim::SciMLBase.EnsembleSolution, success_rate)
        mean_res_time = mean([sol.t[1] for sol in sim])
        mean_trans_time = mean([(sol.t[end] - sol.t[1]) for sol in sim])

        return new{typeof(mean_res_time)}(
            success_rate, mean_res_time, mean_trans_time, mean_res_time / mean_trans_time
        )
    end
end;

"""
$(TYPEDEF)

Ensemble of transition paths between two points in a state space.

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

"""
struct TransitionEnsemble{SSS,T,ES}
    "paths sampled from the transition process"
    paths::Vector{SSS}
    "coresponsing times of the paths"
    times::Vector{Vector{T}}
    "statistics of the ensemble"
    stats::TransitionStatistics{T}
    "original ensemble solution of the SciML Ensemble problem"
    sciml_ensemble::ES

    function TransitionEnsemble(sim::SciMLBase.EnsembleSolution, success_rate)
        stats = TransitionStatistics(sim, success_rate)

        samples = [StateSpaceSet(sol.u) for sol in sim]
        times = [sol.t for sol in sim]

        return new{eltype(samples),eltype(eltype(times)),typeof(sim)}(
            samples, times, stats, sim
        )
    end
end;

function prettyprint(te::TransitionEnsemble)
    ts = te.stats
    return "Transition path ensemble of $(length(te.times)) samples
           - sampling success rate:  $(round(ts.success_rate, sigdigits=3))
           - mean residence time:    $(round(ts.residence_time, sigdigits=3))
           - mean transition time:   $(round(ts.transition_time, sigdigits=3))
           - rareness:               $(round(ts.rareness, sigdigits=3))"
end

Base.show(io::IO, te::TransitionEnsemble) = print(io, prettyprint(te))
