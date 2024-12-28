mutable struct MinimumActionPath{D,T,V,Phis,Ahis,Lambda,PV,GPV}
    path::StateSpaceSet{D,T,V}
    action::T
    path_history::Phis
    action_history::Ahis
    λ::Lambda
    generalized_momentum::GPV
    path_velocity::PV

    function MinimumActionPath(
        path::StateSpaceSet{D,T,V},
        action;
        path_history=nothing,
        action_history=nothing,
        λ=nothing,
        generalized_momentum=nothing,
        path_velocity=nothing,
    ) where {D,T,V}
        return new{
            D,
            T,
            V,
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

function prettyprint(mlp::MinimumActionPath{D}) where {D}
    return "Minimum action Path of length $(length(mlp.path)) in $D dimensions"
end

Base.show(io::IO, mlp::MinimumActionPath) = print(io, prettyprint(mlp))
