"""
    FreidlinWentzellHamiltonian{IIP, D, ...}

Freidlin-Wentzell Hamiltonian for the small-noise SDE
``\\mathrm{d}X_t = b(X_t)\\,\\mathrm{d}t + \\sigma\\,\\Sigma(X_t)\\,\\mathrm{d}W_t``,
with diffusion tensor ``a(x) = \\Sigma(x)\\Sigma(x)^\\top``:

```math
H(x, p) \\;=\\; \\langle b(x),\\, p\\rangle \\;+\\; \\tfrac{1}{2}\\,\\langle p,\\, a(x)\\,p\\rangle.
```

`p` is the conjugate momentum (Legendre dual of ``\\dot x``); Hamilton's equations are
``\\dot x = b(x) + a(x)\\,p`` and ``\\dot p = -\\partial_x b^\\top p - \\tfrac{1}{2}\\langle p,
\\partial_x a\\,p\\rangle``. Freidlin-Wentzell instantons (action minimizers between
invariant sets) live on the zero-energy shell ``H \\equiv 0``; the simplified geometric
MAM ([grafke_long_2017](@cite)) minimizes the action directly in `(x, p)` space.

The stored `a(x)` is trace-normalized (see [`_trace_normalized_a`](@ref)) so that the
action is invariant to the overall scale of `noise_strength`.

## Fields
* `H_x`, `H_p`: callables `(x, p) -> ∂_x H`, `(x, p) -> ∂_p H`, returning `D × N` matrices.
* `a`: trace-normalized diffusion callable; a `Base.Returns` for constant noise.
* `x_ref`: reference state used by `Base.show`; `nothing` if not supplied.

## Constructors
* `FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)`: builds `H_x`, `H_p`,
  `a` from `ds`. Central finite differences are used for ``\\partial_x a`` when `a` is
  state-dependent. Validation and diagonal-vs-coupled classification happen at cache
  build ([`build_sgmam_cache`](@ref)); rank-deficient `a` is rejected there.
* `FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p; a = Returns(Diagonal(ones(D))), x_ref = nothing)`:
  for hand-rolled Hamiltonians; you are responsible for matching the convention above.

See [freidlin_random_1998](@cite) for the underlying theory.
"""
struct FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, A, R}
    H_x::Hx
    H_p::Hp
    a::A
    x_ref::R
end

is_constant(::Base.Returns) = Val(true)
is_constant(_) = Val(false)

function FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    ds isa CoupledSDEs && proper_FW_system(ds)
    a = _trace_normalized_a(ds)
    f = dynamic_rule(ds)
    jac = jacobian(ds)
    ps = current_parameters(ds)
    H_p = _make_H_p(a, f, ps)
    H_x = _make_H_x(a, jac, ps)
    x_ref = collect(current_state(ds))
    return FreidlinWentzellHamiltonian{
        isinplace(ds), D, typeof(H_x), typeof(H_p), typeof(a), typeof(x_ref),
    }(H_x, H_p, a, x_ref)
end

function FreidlinWentzellHamiltonian{IIP, D}(
        H_x, H_p;
        a = Returns(LinearAlgebra.Diagonal(ones(Float64, D))),
        x_ref = nothing,
    ) where {IIP, D}
    return FreidlinWentzellHamiltonian{IIP, D, typeof(H_x), typeof(H_p), typeof(a), typeof(x_ref)}(
        H_x, H_p, a, x_ref,
    )
end

_make_H_p(a, f, ps) = function H_p(x, p)
    Hp = similar(x)
    @inbounds for i in axes(x, 2)
        @views Hp[:, i] .= a(x[:, i]) * p[:, i] .+ f(x[:, i], ps, 0.0)
    end
    return Hp
end

_make_H_x(a, jac, ps) = function H_x(x, p)
    Hx = similar(x)
    e = zeros(eltype(x), size(x, 1))
    @inbounds for i in axes(x, 2)
        @views begin
            LinearAlgebra.mul!(Hx[:, i], jac(x[:, i], ps, 0.0)', p[:, i])
            _add_da_term!(is_constant(a), Hx[:, i], a, x[:, i], p[:, i], e)
        end
    end
    return Hx
end

_add_da_term!(::Val{true}, _, _, _, _, _) = nothing

function _add_da_term!(::Val{false}, out, a, x, p, e)
    h = _fd_step(eltype(x))
    @inbounds for l in eachindex(out)
        e[l] = h
        dla = (a(x .+ e) .- a(x .- e)) ./ (2 * h)
        out[l] += 0.5 * dot(p, dla, p)
        e[l] = 0
    end
    return nothing
end

_a_shape_label(::LinearAlgebra.Diagonal) = "Diagonal a"
_a_shape_label(::AbstractMatrix) = "a"

_a_const_label(::Val{true}) = "constant"
_a_const_label(::Val{false}) = "state-dependent"

function _show_a_label(a, x_ref::AbstractVector)
    sample = a(x_ref)
    return "$(_a_const_label(is_constant(a))) $(_a_shape_label(sample))"
end

_show_a_label(a, ::Nothing) = "$(_a_const_label(is_constant(a))) a"

function Base.show(io::IO, sys::FreidlinWentzellHamiltonian{IIP, D}) where {IIP, D}
    iip = IIP ? "in-place" : "out-of-place"
    label = _show_a_label(sys.a, sys.x_ref)
    return print(
        io,
        "Freidlin-Wentzell Hamiltonian on $D-dimensional state space ($label) with $iip H_x and H_p",
    )
end
