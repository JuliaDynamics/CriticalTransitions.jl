"""
    NoiseShape

Abstract supertype for the diffusion-tensor shape tags used by the Freidlin-Wentzell
machinery. The tag is a singleton type carried as a type parameter so that the per-shape
algorithm branches dispatch at compile time.

Singletons:

* `AdditiveNoise`: `a(x)` is constant (state-independent) and `isdiag(a)` is true.
* `DiagonalNoise`: `a(x)` is diagonal at every probe point but varies with `x`.
* `GeneralNoise`: `a(x)` has off-diagonal entries (constant or state-dependent).

Rank-deficient noise (a singular diffusion tensor) is not supported in this release.
"""
abstract type NoiseShape end

"""
    AdditiveNoise <: NoiseShape

Singleton tag for state-independent diagonal diffusion: `a(x) ≡ const` and `isdiag(a)`.
"""
struct AdditiveNoise <: NoiseShape end

"""
    DiagonalNoise <: NoiseShape

Singleton tag for state-dependent diagonal diffusion: `a(x)` varies with `x` but is
diagonal at every probe point.
"""
struct DiagonalNoise <: NoiseShape end

"""
    GeneralNoise <: NoiseShape

Singleton tag for diffusion with off-diagonal entries, either constant or state-dependent.
"""
struct GeneralNoise <: NoiseShape end

"""
    _is_rank_deficient(M)

Return `true` if the matrix `M` is numerically singular. Uses `cond(M)` against
`1/sqrt(eps(eltype(M)))` as the tolerance.
"""
function _is_rank_deficient(M::AbstractMatrix)
    T = eltype(M)
    tol = 1 / sqrt(eps(real(T)))
    return LinearAlgebra.cond(Matrix(M)) > tol
end

"""
    _probe_points(u₀)

Return `2D + 1` probe points around `u₀` (where `D = length(u₀)`): the center point,
plus `u₀ ± h·eᵢ` for each canonical direction `eᵢ`, with `h = max(√eps, 1e-6)`. Used by
`_classify_noise_shape` to sample `a(x)` for state-dependence and rank-deficiency
detection.
"""
function _probe_points(u₀::AbstractVector)
    D = length(u₀)
    T = eltype(u₀)
    h = max(sqrt(eps(real(T))), T(1.0e-6))
    probes = Vector{Vector{T}}(undef, 2D + 1)
    probes[1] = collect(u₀)
    @inbounds for i in 1:D
        plus = collect(u₀)
        minus = collect(u₀)
        plus[i] += h
        minus[i] -= h
        probes[2i] = plus
        probes[2i + 1] = minus
    end
    return probes
end

"""
    _classify_noise_shape(ds::CoupledSDEs)

Sample `a(x) = σ(x)σ(x)ᵀ` at `_probe_points(current_state(ds))` and return the
appropriate `NoiseShape` singleton. Throws `ArgumentError` if the noise is not
autonomous or if `a` is rank-deficient at any probe.

Uses `diffusion_function(ds)` directly because `covariance_matrix(ds)` returns
`nothing` whenever `:invertible == false` (including for rank-deficient additive
noise, which we want to detect cleanly here).
"""
function _classify_noise_shape(ds::CoupledSDEs)
    ds.noise_type[:autonomous] ||
        throw(ArgumentError("non-autonomous noise not supported"))

    σ_fn = diffusion_function(ds)
    u₀ = current_state(ds)
    ps = current_parameters(ds)

    a_of(u) = let σx = σ_fn(u, ps, 0.0)
        σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
        σ_mat * σ_mat'
    end

    if ds.noise_type[:additive]
        a = a_of(u₀)
        _is_rank_deficient(a) && throw(
            ArgumentError(
                "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p).",
            ),
        )
        return LinearAlgebra.isdiag(a) ? AdditiveNoise() : GeneralNoise()
    else
        probes = _probe_points(u₀)
        as = map(a_of, probes)
        any(_is_rank_deficient, as) && throw(
            ArgumentError(
                "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p).",
            ),
        )
        return all(LinearAlgebra.isdiag, as) ? DiagonalNoise() : GeneralNoise()
    end
end

# `CoupledODEs` has no noise; treat it as identity additive noise for the FW Hamiltonian.
_classify_noise_shape(::CoupledODEs) = AdditiveNoise()

"""
    _classify_user_a(a_callable, D)

Classify a user-supplied `a(x)` callable. Probes `a` around `zeros(D)`; detects
constant-vs-state-dependent structurally (the callable is a `Base.Returns`) rather
than by numerical comparison, which avoids false positives on barely-state-dependent
functions.
"""
function _classify_user_a(a_callable, D::Int)
    probes = _probe_points(zeros(Float64, D))
    as = map(a_callable, probes)
    if any(_is_rank_deficient, as)
        throw(
            ArgumentError(
                "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p).",
            ),
        )
    end
    if a_callable isa Base.Returns
        return LinearAlgebra.isdiag(as[1]) ? AdditiveNoise() : GeneralNoise()
    end
    return all(LinearAlgebra.isdiag, as) ? DiagonalNoise() : GeneralNoise()
end
