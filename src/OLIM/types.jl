
struct MyCurve
    x::Float64
    y::Float64
    next::Union{MyCurve, Nothing}
    prev::Union{MyCurve, Nothing}
end

struct MySol
    g::Float64
    c::Char
end

struct SolInfo
    type::Char # solution type: 1 = 1ptupdate, 2 = 2ptupdate, 0 = initialization, 'n' = never reached
    ind0::Int
    ind1::Int
    s::Float64
end
