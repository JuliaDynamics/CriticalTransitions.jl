# LimitCycleFrame: moving frame and normal-plane coefficients along a stable limit cycle.

using DynamicalSystemsBase: TangentDynamicalSystem, current_deviations
using LinearAlgebra: I, eigen, norm
using StateSpaceSets: StateSpaceSet

"""
    _state_transition_matrices(sys::CoupledSDEs, points::StateSpaceSet, period::Real)

Return `Φ::Array{T,3}` of size `d × d × Nτ` with `Φ[:, :, k]` equal to the state
transition matrix `Φ(τ_k, 0)` of the variational equation `Φ̇ = ∂ₓb(γ(τ))·Φ` along the
periodic orbit, with `Φ(0,0) = I`.
"""
function _state_transition_matrices(sys::CoupledSDEs, points::StateSpaceSet, period::Real)
    d = length(points[1])
    Nτ = length(points)
    Δτ = period / Nτ
    Tf = eltype(points[1])
    Φs = zeros(Tf, d, d, Nτ)
    Φs[:, :, 1] = Matrix{Tf}(I, d, d)
    tands = TangentDynamicalSystem(CoupledODEs(sys); u0 = collect(points[1]))
    for k in 2:Nτ
        DynamicalSystemsBase.step!(tands, Δτ, true)
        Φs[:, :, k] = current_deviations(tands)
    end
    return Φs
end

"""
    LimitCycleFrame{D, T}

Moving frame and normal-plane coefficient matrices along a stable limit cycle Γ.
"""
struct LimitCycleFrame{D, T}
    γ::StateSpaceSet{D, T}
    period::T
    E::Array{T, 3}
    Ẽ::Array{T, 3}
    M̃::Array{T, 3}
    Ã::Array{T, 3}
end

function LimitCycleFrame(
        points::StateSpaceSet{D, T}, period::Real, sys::CoupledSDEs;
        tol_periodic::Real = 1.0e-3,
    ) where {D, T}
    period_t = T(period)
    Φs = _state_transition_matrices(sys, points, period_t)
    E = _build_frame(sys, points, period_t, Φs; tol_periodic)
    Ẽ = E[:, 2:D, :]
    M̃, Ã = _assemble_M_A(sys, points, period_t, E, Ẽ)
    return LimitCycleFrame{D, T}(points, period_t, E, Ẽ, M̃, Ã)
end

function _build_frame(sys, points, period, Φs; tol_periodic)
    d = length(points[1])
    Nτ = length(points)
    E = zeros(d, d, Nτ)
    for k in 1:Nτ
        v = drift(sys, collect(points[k]))
        E[:, 1, k] .= v ./ norm(v)
    end
    tands = TangentDynamicalSystem(CoupledODEs(sys); u0 = collect(points[1]))
    DynamicalSystemsBase.step!(tands, period, true)
    Φ_M = Matrix(current_deviations(tands))
    F = eigen(Matrix(Φ_M'))
    λs = F.values
    V = F.vectors
    trivial = argmin(abs.(λs .- 1))
    nontrivial = [i for i in 1:d if i != trivial]
    e_initial = _real_basis(V[:, nontrivial], λs[nontrivial])
    # Sign-continue each transported vector against its predecessor; otherwise successive
    # `\` solves can produce arbitrary sign flips that look like anti-periodicity.
    for j in 1:(d - 1)
        v0 = e_initial[:, j] / norm(e_initial[:, j])
        E[:, j + 1, 1] .= v0
        prev = v0
        for k in 2:Nτ
            v = transpose(Φs[:, :, k]) \ e_initial[:, j]
            v ./= norm(v)
            if dot(v, prev) < 0
                v .= -v
            end
            E[:, j + 1, k] .= v
            prev = v
        end
    end
    # Per-vector orientability at τ = T. `rel = min(‖v_T - v_0‖, ‖v_T + v_0‖)` collapses
    # sign, so both anti-periodic (`v_T ≈ -v_0`) and non-orientable cases trip one threshold.
    for j in 1:(d - 1)
        v_T = transpose(Φ_M) \ e_initial[:, j]
        v_T ./= norm(v_T)
        v_0 = E[:, j + 1, 1]
        rel = min(norm(v_T - v_0), norm(v_T + v_0))
        if rel > tol_periodic
            throw(
                ArgumentError(
                    "LimitCycleFrame: vector j=$(j) is not T-periodic after transport " *
                        "(residual = $(rel)). Anti-periodic or non-orientable frames are not supported.",
                )
            )
        end
    end
    return E
end

# Real (d-1)-column basis from complex eigenvectors; conjugate pairs → 2D real blocks.
function _real_basis(Vs::AbstractMatrix, λs::AbstractVector)
    d = size(Vs, 1)
    m = length(λs)
    out = zeros(d, m)
    used = falses(m)
    col = 1
    for i in 1:m
        used[i] && continue
        if isreal(λs[i])
            out[:, col] = real.(Vs[:, i])
            used[i] = true
            col += 1
        else
            j = findfirst(k -> !used[k] && k != i && isapprox(λs[k], conj(λs[i])), 1:m)
            j === nothing && error("Unpaired complex eigenvalue in monodromy basis.")
            out[:, col] = real.(Vs[:, i])
            out[:, col + 1] = imag.(Vs[:, i])
            used[i] = true
            used[j] = true
            col += 2
        end
    end
    return out
end

function _assemble_M_A(sys, points, period, E, Ẽ)
    d = length(points[1])
    Nτ = length(points)
    Δτ = period / Nτ
    a_fn = _trace_normalized_a(sys)
    drift_fn = let p = current_parameters(sys), f = dynamic_rule(sys)
        x -> f(x, p, 0.0)
    end
    M̃ = zeros(d - 1, d - 1, Nτ)
    Ã = zeros(d - 1, d - 1, Nτ)
    Erecip = zeros(d, d, Nτ)
    for k in 1:Nτ
        Erecip[:, :, k] = inv(E[:, :, k])
    end
    Ė = similar(E)
    for k in 1:Nτ
        kp = mod1(k + 1, Nτ)
        km = mod1(k - 1, Nτ)
        Ė[:, :, k] = (E[:, :, kp] - E[:, :, km]) / (2 * Δτ)
    end
    for k in 1:Nτ
        xk = collect(points[k])
        J = ForwardDiff.jacobian(drift_fn, xk)
        a_k = a_fn(xk)
        for i in 1:(d - 1), j in 1:(d - 1)
            J̃ij = dot(Erecip[i + 1, :, k], J * Ẽ[:, j, k])
            Ωij = dot(Erecip[i + 1, :, k], Ė[:, j + 1, k])
            M̃[i, j, k] = J̃ij - Ωij
            Ã[i, j, k] = dot(Erecip[i + 1, :, k], a_k * Erecip[j + 1, :, k])
        end
    end
    return M̃, Ã
end
