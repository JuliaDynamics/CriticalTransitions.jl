struct TransitionStatistics{T}
    success_rate::T
    residence_time::T
    transition_time::T
    rareness::T

    function TransitionStatistics(sim, success_rate)
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
struct TransitionEnsemble{SSS, T, ES}
    paths::Vector{SSS}
    times::Vector{Vector{T}}
    stats::TransitionStatistics{T}
    sciml_ensemble::ES

    function TransitionEnsemble(sim, success_rate)
        stats = TransitionStatistics(sim, success_rate)

        samples = [StateSpaceSet(sol.u) for sol in sim]
        times = [sol.t for sol in sim]

        return new{eltype(samples),eltype(eltype(times)),typeof(sim)}(samples, times, stats, sim)
    end
end;

function prettyprint(te::TransitionEnsemble)
    ts = te.stats
    return "Transition path ensemble of $(length(te.times)) samples
           - sampling success rate:      $(round(ts.success_rate, digits=3))
           - mean residence time:        $(round(ts.residence_time, digits=3))
           - mean transition time:       $(round(ts.transition_time, digits=3))
           - rareness: $(round(ts.rareness, digits=1))"
end

Base.show(io::IO, te::TransitionEnsemble) = print(io, prettyprint(te))
