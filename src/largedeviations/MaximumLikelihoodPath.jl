mutable struct MaximumLikelihoodPath{T,Phis,Ahis,Lambda,PV,GPV}
    path::Matrix{T}
    action::T
    path_history::Phis
    action_history::Ahis
    λ::Lambda
    generalized_momentum::GPV
    path_velocity::PV

    function MaximumLikelihoodPath(
        path::Matrix{T},
        action;
        path_history=nothing,
        action_history=nothing,
        λ=nothing,
        generalized_momentum=nothing,
        path_velocity=nothing,
    ) where {T}
        return new{
            T,
            typeof(path_history),
            typeof(action_history),
            typeof(λ),
            typeof(generalized_momentum),
            typeof(path_velocity),
        }(
            path,
            action,
            path_history,
            action_history,
            λ,
            generalized_momentum,
            path_velocity,
        )
    end
end
