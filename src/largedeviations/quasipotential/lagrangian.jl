"""
    itp_root(f, a, b; k1, k2, n0, atol, maxiter) -> (root::T, converged::Bool)

Interpolate-Truncate-Project bracketed scalar root-find (Oliveira and Takahashi 2020,
ACM TOMS 47:1). Returns `(root, true)` if a sign change is bracketed; returns the
endpoint with smaller `|f|` and `false` when `f(a) * f(b) > 0`, signalling no
interior root.
"""
@inline function itp_root(
        f, a::T, b::T;
        k1::T = T(0.1), k2::T = T(2),
        n0::Int = 1,
        atol::T = eps(T)^(T(3) / T(4)),
        maxiter::Int = 64,
    ) where {T <: AbstractFloat}
    fa = f(a); fb = f(b)
    iszero(fa) && return (a, true)
    iszero(fb) && return (b, true)
    if fa * fb > 0
        return (abs(fa) <= abs(fb) ? a : b, false)
    end
    ε = T(2) * atol
    n12 = ceil(Int, log2((b - a) / ε))
    nmax = n12 + n0
    aₖ, bₖ = a, b
    faₖ, fbₖ = fa, fb
    for j in 0:(maxiter - 1)
        xf = (faₖ * bₖ - fbₖ * aₖ) / (faₖ - fbₖ)
        xh = (aₖ + bₖ) / 2
        σ = sign(xh - xf)
        δ = k1 * (bₖ - aₖ)^k2
        xt = δ <= abs(xh - xf) ? xf + σ * δ : xh
        r = ε * T(2)^(T(nmax) - T(j)) - (bₖ - aₖ) / 2
        xp = abs(xt - xh) <= r ? xt : xh - σ * r
        fp = f(xp)
        if iszero(fp)
            return (xp, true)
        elseif fp * faₖ < 0
            bₖ, fbₖ = xp, fp
        else
            aₖ, faₖ = xp, fp
        end
        bₖ - aₖ < ε && break
    end
    return ((aₖ + bₖ) / 2, true)
end

"""
Internal wrapper bundling drift, trace-normalized inverse metric, and the
near-zero-drift threshold for `L_g(x, v) = |v|_Q |b|_Q - <v, b>_Q`.
`Q_inv` is either an `SMatrix{D, D, T, D*D}` (additive noise) or a callable
`x -> SMatrix{D, D, T, D*D}` (multiplicative noise).
"""
struct _GeometricLagrangian{D, T, B, A}
    b::B
    Q_inv::A
    eps_b::T
end

_GeometricLagrangian{D, T}(b::B, Q_inv::A, eps_b::T) where {D, T, B, A} =
    _GeometricLagrangian{D, T, B, A}(b, Q_inv, eps_b)

@inline _Q_inv_at(M::SMatrix, _) = M
@inline _Q_inv_at(f::F, x) where {F} = f(x)

@inline function (L::_GeometricLagrangian{D, T})(
        x::SVector{D, S}, v::SVector{D, S}
    ) where {D, T, S}
    bx = L.b(x)
    Qinv = _Q_inv_at(L.Q_inv, x)
    bnQ2 = dot(bx, Qinv * bx)
    if bnQ2 < S(L.eps_b) * S(L.eps_b)
        return S(0.5) * dot(v, Qinv * v) - dot(v, Qinv * bx)
    end
    return sqrt(dot(v, Qinv * v)) * sqrt(bnQ2) - dot(v, Qinv * bx)
end

struct _DriftOOP{D, F, P}
    raw::F
    p::P
end
@inline function (b::_DriftOOP{D})(x::SVector{D, S}) where {D, S}
    return SVector{D, S}(b.raw(x, b.p, zero(S)))
end
_make_b_oop(raw, p, ::Val{D}) where {D} = _DriftOOP{D, typeof(raw), typeof(p)}(raw, p)

struct _DriftIIP{D, F, P}
    raw::F
    p::P
end
@inline function (b::_DriftIIP{D})(x::SVector{D, S}) where {D, S}
    buf = MVector{D, S}(undef)
    b.raw(buf, x, b.p, zero(S))
    return SVector{D, S}(buf)
end
_make_b_iip(raw, p, ::Val{D}) where {D} = _DriftIIP{D, typeof(raw), typeof(p)}(raw, p)

struct _QInvDynamic{D, T, F}
    a_fn::F
end
@inline function (q::_QInvDynamic{D, T})(x::SVector{D}) where {D, T}
    return inv(SMatrix{D, D, T}(q.a_fn(x)))
end
_make_Q_inv(a_fn, ::Val{D}, ::Type{T}) where {D, T} =
    _QInvDynamic{D, T, typeof(a_fn)}(a_fn)

"""
Build a `_GeometricLagrangian` from `sys`. Drift via `dynamic_rule(sys)`,
wrapped to return `SVector{D, T}` regardless of `IIP`. `Q_inv` via
`_trace_normalized_a(sys)`, inverted to an `SMatrix` (additive) or wrapped
as `x -> inv(SMatrix(...))` (multiplicative).
"""
function _geometric_lagrangian(
        sys::CoupledSDEs{IIP, D}, ::Type{T};
        eps_b::T = T(1.0e-10),
    ) where {IIP, D, T}
    raw = dynamic_rule(sys); p = current_parameters(sys)
    b = if IIP
        _make_b_iip(raw, p, Val(D))
    else
        _make_b_oop(raw, p, Val(D))
    end
    a_fn = _trace_normalized_a(sys)
    Q_inv = if a_fn isa Base.Returns
        inv(SMatrix{D, D, T}(a_fn(current_state(sys))))
    else
        _make_Q_inv(a_fn, Val(D), T)
    end
    return _GeometricLagrangian{D, T}(b, Q_inv, eps_b)
end

@inline function _Lg_at(
        L::_GeometricLagrangian{D, T},
        x::SVector{D, S}, Qinv_x::SMatrix{D, D},
        Qv::SVector{D}, vQv::S, sqrt_vQv::S,
    ) where {D, T, S}
    bx = L.b(x)
    bnQ2 = dot(bx, Qinv_x * bx)
    if bnQ2 < S(L.eps_b) * S(L.eps_b)
        return S(0.5) * vQv - dot(Qv, bx)
    end
    return sqrt_vQv * sqrt(bnQ2) - dot(Qv, bx)
end

# Additive noise: Qinv is a constant SMatrix; hoist Qv, vQv, sqrt(vQv) out
# of the three Simpson evaluations.
@inline function _line_integral(
        L::_GeometricLagrangian{D, T, B, A},
        y::SVector{D, S}, v::SVector{D, S},
    ) where {D, T, B, A <: SMatrix, S}
    half = S(0.5)
    Qinv = L.Q_inv
    Qv = Qinv * v
    vQv = dot(v, Qv)
    sv = sqrt(vQv)
    L0 = _Lg_at(L, y, Qinv, Qv, vQv, sv)
    Lh = _Lg_at(L, y + half * v, Qinv, Qv, vQv, sv)
    L1 = _Lg_at(L, y + v, Qinv, Qv, vQv, sv)
    return (L0 + S(4) * Lh + L1) / S(6)
end

# Multiplicative noise fallback (Qinv depends on x).
@inline function _line_integral(
        L::_GeometricLagrangian{D, T}, y::SVector{D, S}, v::SVector{D, S},
    ) where {D, T, S}
    half = S(0.5)
    return (L(y, v) + S(4) * L(y + half * v, v) + L(y + v, v)) / S(6)
end

@inline function _hermite_U(U0::T, U1::T, m0::T, m1::T, λ::T) where {T}
    s0 = isnan(m0) ? (U1 - U0) : m0
    s1 = isnan(m1) ? (U1 - U0) : m1
    h00 = T(2) * λ^3 - T(3) * λ^2 + one(T)
    h10 = λ^3 - T(2) * λ^2 + λ
    h01 = -T(2) * λ^3 + T(3) * λ^2
    h11 = λ^3 - λ^2
    return h00 * U0 + h10 * s0 + h01 * U1 + h11 * s1
end
