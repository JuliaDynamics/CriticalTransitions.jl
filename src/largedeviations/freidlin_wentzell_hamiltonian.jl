"""
A structure representing an extended phase space system where the dissipative vector
field is embedded in a doubled dimensional phase space. Given the original phase space
coordinates `x` of a vector field `b(x)`, we define the canonical momenta `p` so that
the new phase space coordinates are `(x, p)`. The dynamics in this extended phase
space are governed by the Freidlin/Wentzell Hamiltonian:

``H(x, p) = ⟨b(x), p⟩ + (1/2) ⟨p, a(x)·p⟩,``

where `a(x) = σ(x)σ(x)ᵀ` is the diffusion tensor. The Hamiltonian is assumed to be
quadratic in `p`; geodesics on the zero-energy shell `H ≡ 0` are Freidlin-Wentzell
instantons.

The struct stores `H_x` and `H_p` (partial derivatives of `H` with respect to `x` and
`p`), a callable `a` for the diffusion tensor, and a [`NoiseShape`](@ref) type
parameter `NS` that encodes how the inner algorithm loops dispatch:
* `AdditiveNoise`: constant diagonal `a`.
* `DiagonalNoise`: state-dependent diagonal `a(x)`.
* `GeneralNoise`: full-matrix `a(x)`, constant or state-dependent.

# Constructors
* `FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)`: builds `H_x`,
  `H_p`, and the trace-normalized `a` from the system's drift and diffusion. The
  resulting `H_x` differentiates `a(x)` by finite differences when the noise is
  state-dependent.
* `FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p; a = Returns(I))`: user-supplied
  analytic Hamilton's functions. Pass an `a` callable to indicate a non-identity
  diffusion tensor; the constructor classifies its structure into a `NoiseShape`.
"""
struct FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS <: NoiseShape}
    H_x::Hx
    H_p::Hp
    a::AF
end

function FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    if ds isa CoupledSDEs
        proper_FW_system(ds)
    end
    NS = _classify_noise_shape(ds)

    a_callable = _trace_normalized_a(ds)
    ps = current_parameters(ds)

    is_constant_a = (ds isa CoupledODEs) || (ds isa CoupledSDEs && ds.noise_type[:additive])

    f = dynamic_rule(ds)
    jac = jacobian(ds)

    H_x = _build_H_x(NS, is_constant_a, jac, a_callable, ps)

    function H_p(x, p)
        Hp = similar(x)
        @views @inbounds for idx in 1:size(x, 2)
            xi = x[:, idx]; pi_v = p[:, idx]
            Hp[:, idx] .= a_callable(xi) * pi_v .+ f(xi, ps, 0.0)
        end
        return Hp
    end

    return FreidlinWentzellHamiltonian{
        isinplace(ds), D, typeof(H_x), typeof(H_p), typeof(a_callable), typeof(NS),
    }(H_x, H_p, a_callable)
end

function _build_H_x(NS, is_constant_a, jac, a_callable, ps)
    if is_constant_a
        return _hx_drift_only(jac, ps)
    elseif NS isa DiagonalNoise
        return _hx_diagonal(jac, ps, a_callable)
    else
        return _hx_general(jac, ps, a_callable)
    end
end

function _hx_drift_only(jac, ps)
    return function H_x_drift(x, p)
        Hx = similar(x)
        Nt = size(x, 2); Dn = size(x, 1)
        @views @inbounds for idx in 1:Nt
            xi = x[:, idx]; pi_v = p[:, idx]
            jax = jac(xi, ps, 0.0)
            for idc in 1:Dn
                Hx[idc, idx] = dot(jax[:, idc], pi_v)
            end
        end
        return Hx
    end
end

function _hx_diagonal(jac, ps, a_callable)
    return function H_x_diag(x, p)
        Hx = similar(x)
        Nt = size(x, 2); Dn = size(x, 1)
        h_fd = _fd_step(eltype(x))
        e = zeros(eltype(x), Dn)
        @views @inbounds for idx in 1:Nt
            xi = x[:, idx]; pi_v = p[:, idx]
            jax = jac(xi, ps, 0.0)
            for idc in 1:Dn
                Hx[idc, idx] = dot(jax[:, idc], pi_v)
            end
            for l in 1:Dn
                dla_diag = LinearAlgebra.diag(_da_dx_l(a_callable, xi, l, h_fd, e))
                acc = zero(eltype(x))
                for k in 1:Dn
                    acc += pi_v[k]^2 * dla_diag[k]
                end
                Hx[l, idx] += 0.5 * acc
            end
        end
        return Hx
    end
end

function _hx_general(jac, ps, a_callable)
    return function H_x_gen(x, p)
        Hx = similar(x)
        Nt = size(x, 2); Dn = size(x, 1)
        h_fd = _fd_step(eltype(x))
        e = zeros(eltype(x), Dn)
        @views @inbounds for idx in 1:Nt
            xi = x[:, idx]; pi_v = p[:, idx]
            jax = jac(xi, ps, 0.0)
            for idc in 1:Dn
                Hx[idc, idx] = dot(jax[:, idc], pi_v)
            end
            for l in 1:Dn
                dla = _da_dx_l(a_callable, xi, l, h_fd, e)
                acc = zero(eltype(x))
                for k1 in 1:Dn, k2 in 1:Dn
                    acc += pi_v[k1] * dla[k1, k2] * pi_v[k2]
                end
                Hx[l, idx] += 0.5 * acc
            end
        end
        return Hx
    end
end

function FreidlinWentzellHamiltonian{IIP, D}(
        H_x::Function, H_p::Function;
        a = Returns(LinearAlgebra.Diagonal(ones(Float64, D))),
    ) where {IIP, D}
    NS = _classify_user_a(a, D)
    return FreidlinWentzellHamiltonian{IIP, D, typeof(H_x), typeof(H_p), typeof(a), typeof(NS)}(
        H_x, H_p, a,
    )
end

function prettyprint(mlp::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS}) where {IIP, D, Hx, Hp, AF, NS}
    return "Freidlin-Wentzell Hamiltonian on $D-dimensional state space ($(NS)) with $(IIP ? "in-place" : "out-of-place") H_x and H_p"
end

Base.show(io::IO, mlp::FreidlinWentzellHamiltonian) = print(io, prettyprint(mlp))
