"""
$(TYPEDEF)

The minimum action path between two points in a D-dimensional phase space.

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

"""
struct MinimumActionPath{D,T<:Real,V,Phis,Ahis,Lambda,PV,GPV}
    """The path matrix."""
    path::StateSpaceSet{D,T,V}
    """The action value associated to the path."""
    action::T
    """The history of action of the paths in the optimisation algorithm (optional)."""
    path_history::Phis
    """The history of action of the paths in the optimisation algorithm (optional)."""
    action_history::Ahis
    """The Lagrange multiplier parameter for the minimum action path (optional)."""
    λ::Lambda
    """The generalized momentum of the phase space variables (optional)."""
    generalized_momentum::GPV
    """The path velocity (optional)."""
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
