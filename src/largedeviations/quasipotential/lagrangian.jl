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

const _OLIM_MULTI_NOISELESS_MSG =
    "OLIM supports at most one noiseless coordinate (|Z| = 1). This diffusion has " *
    "two or more zero coordinates (a sub-Riemannian problem), which is not supported. " *
    _RANK_DEFICIENT_MSG

const _OLIM_ROTATED_NULL_MSG =
    "OLIM requires either an invertible diffusion or a single coordinate-aligned " *
    "noiseless coordinate. This diffusion is singular with a non-coordinate-aligned " *
    "null space, which is not supported. " * _RANK_DEFICIENT_MSG

"""
    _degenerate_split(a0; atol) -> (is_degenerate, z, R_idx)

Classify the trace-normalized diffusion `a0`. Returns `(false, 0, 1:D)` for full rank,
`(true, z, R_idx)` when exactly one coordinate `z` is noiseless (its row and column are
numerically zero) and the remaining block `a0[R, R]` is positive-definite. Throws for two
or more noiseless coordinates or a singular `a0` with a non-coordinate-aligned null space.
"""
function _degenerate_split(a0::AbstractMatrix{T}; atol::T = sqrt(eps(real(T)))) where {T}
    D = LinearAlgebra.checksquare(a0)
    is_zero_rowcol(i) = all(j -> abs(a0[i, j]) <= atol, 1:D) &&
        all(j -> abs(a0[j, i]) <= atol, 1:D)
    Z = findall(is_zero_rowcol, 1:D)
    if isempty(Z)
        if LinearAlgebra.cond(Matrix(a0)) > 1 / sqrt(eps(real(T)))
            throw(ArgumentError(_OLIM_ROTATED_NULL_MSG))
        end
        return (false, 0, SVector{D, Int}(ntuple(identity, D)))
    elseif length(Z) == 1
        z = Z[1]
        Ridx = SVector{D - 1, Int}(ntuple(k -> k < z ? k : k + 1, D - 1))
        A_RR = Matrix(a0[Ridx, Ridx])
        LinearAlgebra.isposdef(LinearAlgebra.Symmetric(A_RR)) ||
            throw(ArgumentError(_OLIM_ROTATED_NULL_MSG))
        return (true, z, Ridx)
    else
        throw(ArgumentError(_OLIM_MULTI_NOISELESS_MSG))
    end
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

@inline _onehot_diag(::Val{D}, z::Int, ::Type{T}) where {D, T} =
    SMatrix{D, D, T}(LinearAlgebra.Diagonal(SVector{D, T}(ntuple(i -> T(i == z), D))))

# Multiplicative-noise inverse metric with the noiseless coordinate `z` regularized by
# `delta`, so the otherwise-singular trace-normalized diffusion is invertible. The
# `delta * e_z e_z'` regularizer is constant, so it is precomputed once into `reg`.
struct _QInvRegDynamic{D, T, F, L}
    a_fn::F
    reg::SMatrix{D, D, T, L}
end
@inline function (q::_QInvRegDynamic{D, T})(x::SVector{D}) where {D, T}
    a = SMatrix{D, D, T}(q.a_fn(x))
    return inv(a + q.reg)
end
_make_Q_inv_reg(a_fn, z::Int, delta::T, ::Val{D}, ::Type{T}) where {D, T} =
    _QInvRegDynamic{D, T, typeof(a_fn), D * D}(a_fn, delta * _onehot_diag(Val(D), z, T))

# Diffusion tensor `a(x)` for the analytic CARE seed: invert the stored `Q_inv`.
@inline _diffusion_at(L::_GeometricLagrangian, x) = inv(_Q_inv_at(L.Q_inv, x))

"""
Build the per-segment Lagrangian from `sys`. Drift via `dynamic_rule(sys)`, wrapped to
return `SVector{D, T}` regardless of `IIP`. The trace-normalized diffusion is classified
by `_degenerate_split`: full rank builds a `_GeometricLagrangian` directly; a single
noiseless coordinate is regularized (`a + regularization * e_z e_z'`, leaving the noisy
block exact) and then builds a `_GeometricLagrangian`; unsupported degeneracy throws.
"""
function _geometric_lagrangian(
        sys::CoupledSDEs{IIP, D}, ::Type{T};
        eps_b::T = T(1.0e-10), regularization::Real = zero(T),
    ) where {IIP, D, T}
    raw = dynamic_rule(sys); p = current_parameters(sys)
    b = if IIP
        _make_b_iip(raw, p, Val(D))
    else
        _make_b_oop(raw, p, Val(D))
    end
    a_fn = _trace_normalized_a(sys)
    a0 = SMatrix{D, D, T}(a_fn(current_state(sys)))
    is_deg, z, _ = _degenerate_split(a0)
    if !is_deg
        Q_inv = a_fn isa Base.Returns ? inv(a0) : _make_Q_inv(a_fn, Val(D), T)
        return _GeometricLagrangian{D, T}(b, Q_inv, eps_b)
    end
    δ = T(regularization)
    δ > zero(T) || throw(
        ArgumentError(
            "rank-1 diffusion needs regularization > 0 to invert the metric; got $δ",
        ),
    )
    Q_inv = a_fn isa Base.Returns ?
        inv(a0 + δ * _onehot_diag(Val(D), z, T)) :
        _make_Q_inv_reg(a_fn, z, δ, Val(D), T)
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

# Precomputed drift data at a cell center: `b = b(c)`, `q = |b|²_Q`, `sq = |b|_Q`.
# Lets the additive line integral reuse the s=0/s=1 Simpson nodes that always land
# on a cell center, skipping a drift eval, a matvec, and a sqrt per cached node.
struct _NodeData{D, T}
    b::SVector{D, T}
    q::T
    sq::T
end

@inline function _Lg_cached(
        eps_b::T, nd::_NodeData{D, T},
        Qv::SVector{D}, vQv::S, sqrt_vQv::S,
    ) where {D, T, S}
    if nd.q < S(eps_b) * S(eps_b)
        return S(0.5) * vQv - dot(Qv, nd.b)
    end
    return sqrt_vQv * nd.sq - dot(Qv, nd.b)
end

# A Simpson node is either cached (`_NodeData`) or evaluated live (`nothing`).
@inline _node_Lg(L::_GeometricLagrangian, x, ::Nothing, Qinv, Qv, vQv, sqrt_vQv) =
    _Lg_at(L, x, Qinv, Qv, vQv, sqrt_vQv)
@inline _node_Lg(L::_GeometricLagrangian, _x, nd::_NodeData, _Qinv, Qv, vQv, sqrt_vQv) =
    _Lg_cached(L.eps_b, nd, Qv, vQv, sqrt_vQv)

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

# Additive noise with cached endpoint nodes. `nd0`/`nd1` are the s=0/s=1 Simpson
# nodes: a `_NodeData` (cell-center lookup) or `nothing` (evaluate live). The
# midpoint is always live. Bitwise-identical to the uncached path.
@inline function _line_integral(
        L::_GeometricLagrangian{D, T, B, A},
        y::SVector{D, S}, v::SVector{D, S}, nd0, nd1,
    ) where {D, T, B, A <: SMatrix, S}
    half = S(0.5)
    Qinv = L.Q_inv
    Qv = Qinv * v
    vQv = dot(v, Qv)
    sv = sqrt(vQv)
    L0 = _node_Lg(L, y, nd0, Qinv, Qv, vQv, sv)
    Lh = _node_Lg(L, y + half * v, nothing, Qinv, Qv, vQv, sv)
    L1 = _node_Lg(L, y + v, nd1, Qinv, Qv, vQv, sv)
    return (L0 + S(4) * Lh + L1) / S(6)
end

# Multiplicative noise has no constant metric to cache against; ignore the nodes.
@inline _line_integral(
    L::_GeometricLagrangian{D, T}, y::SVector{D, S}, v::SVector{D, S}, _nd0, _nd1,
) where {D, T, S} = _line_integral(L, y, v)

@inline function _hermite_U(U0::T, U1::T, m0::T, m1::T, λ::T) where {T}
    s0 = isnan(m0) ? (U1 - U0) : m0
    s1 = isnan(m1) ? (U1 - U0) : m1
    h00 = T(2) * λ^3 - T(3) * λ^2 + one(T)
    h10 = λ^3 - T(2) * λ^2 + λ
    h01 = -T(2) * λ^3 + T(3) * λ^2
    h11 = λ^3 - λ^2
    return h00 * U0 + h10 * s0 + h01 * U1 + h11 * s1
end
