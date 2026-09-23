@doc raw"""
    FreidlinWentzellHamiltonian{IIP, D, ...}

Freidlin-Wentzell Hamiltonian for the small-noise SDE
``\mathrm{d}X_t = b(X_t)\,\mathrm{d}t + \sigma\,\Sigma(X_t)\,\mathrm{d}W_t``,
with diffusion tensor ``a(x) = \Sigma(x)\Sigma(x)^\top``:

```math
H(x, p) \;=\; \langle b(x),\, p\rangle \;+\; \tfrac{1}{2}\,\langle p,\, a(x)\,p\rangle.
```

`p` is the conjugate momentum (Legendre dual of ``\dot x``); Hamilton's equations are
``\dot x = b(x) + a(x)\,p`` and ``\dot p = -\partial_x b^\top p - \tfrac{1}{2}\langle p,
\partial_x a\,p\rangle``. Freidlin-Wentzell instantons (action minimizers between
invariant sets) live on the zero-energy shell ``H \equiv 0``; the simplified geometric
MAM ([grafke_long_2017](@cite)) minimizes the action directly in `(x, p)` space.

The stored `a(x)` is trace-normalized (internal helper `_trace_normalized_a`) so that
the action is invariant to the overall scale of `noise_strength`.

## Fields
* `H_x`, `H_p`: Hamiltonian derivative callables. In-place callables support
  `(out, x, p)`; auto-derived callables also retain the allocating `(x, p)` form for
  compatibility.
* `a`: trace-normalized diffusion callable; a `Base.Returns` for constant noise.
* `x_ref`: reference state used by `Base.show`; `nothing` if not supplied.

## Constructors
* `FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)`: builds in-place
  `H_x`, `H_p`, and `a` from `ds`. Central finite differences are used for
  ``\partial_x a`` when `a` is state-dependent. Validation and diagonal-vs-coupled
  classification happen at cache build (internal helper `build_sgmam_cache`);
  rank-deficient `a` is rejected there.
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

struct _AutoHamiltonianPOOP{A, F, P}
    a::A
    f::F
    ps::P
end

struct _AutoHamiltonianPIIP{A, F, P}
    a::A
    f::F
    ps::P
end

struct _AutoHamiltonianXOOP{A, F}
    a::A
    jac_fn::F
end

struct _AutoHamiltonianXJac{A, J, P}
    a::A
    jac::J
    ps::P
end

function (H::_AutoHamiltonianPOOP)(out, x, p)
    @inbounds for i in axes(x, 2)
        xi = view(x, :, i)
        pi = view(p, :, i)
        oi = view(out, :, i)
        LinearAlgebra.mul!(oi, H.a(xi), pi)
        fi = H.f(xi, H.ps, 0.0)
        for k in eachindex(oi)
            oi[k] += fi[k]
        end
    end
    return out
end

function (H::_AutoHamiltonianPIIP)(out, x, p)
    @inbounds for i in axes(x, 2)
        xi = view(x, :, i)
        pi = view(p, :, i)
        oi = view(out, :, i)
        H.f(oi, xi, H.ps, 0.0)
        ai = H.a(xi)
        for k in eachindex(oi)
            aikp = zero(eltype(out))
            for l in eachindex(pi)
                aikp += ai[k, l] * pi[l]
            end
            oi[k] += aikp
        end
    end
    return out
end

function (H::_AutoHamiltonianXOOP)(out, x, p)
    x_buf = collect(view(x, :, first(axes(x, 2))))
    J_buf = Matrix{eltype(out)}(undef, length(x_buf), length(x_buf))
    jac_cfg = ForwardDiff.JacobianConfig(H.jac_fn, x_buf)
    x_probe = similar(x_buf)
    @inbounds for i in axes(x, 2)
        xi = view(x, :, i)
        pi = view(p, :, i)
        oi = view(out, :, i)
        copyto!(x_buf, xi)
        ForwardDiff.jacobian!(J_buf, H.jac_fn, x_buf, jac_cfg, Val{false}())
        LinearAlgebra.mul!(oi, J_buf', pi)
        _add_da_term!(is_constant(H.a), oi, H.a, xi, pi, x_probe)
    end
    return out
end

function (H::_AutoHamiltonianXJac)(out, x, p)
    x_probe = collect(view(x, :, first(axes(x, 2))))
    @inbounds for i in axes(x, 2)
        xi = view(x, :, i)
        pi = view(p, :, i)
        oi = view(out, :, i)
        LinearAlgebra.mul!(oi, H.jac(xi, H.ps, 0.0)', pi)
        _add_da_term!(is_constant(H.a), oi, H.a, xi, pi, x_probe)
    end
    return out
end

function (H::Union{_AutoHamiltonianPOOP, _AutoHamiltonianPIIP, _AutoHamiltonianXOOP, _AutoHamiltonianXJac})(x, p)
    out = similar(x)
    H(out, x, p)
    return out
end

raw"""
In-place evaluation of ``\partial_p H``. For `IIP=true`, calls `sys.H_p(buf, x, p)`
(user-supplied or auto-derived in-place form). For `IIP=false`, calls `sys.H_p(x, p)`
and copies the returned matrix into `buf`.
"""
@inline _eval_Hp!(buf, sys::FreidlinWentzellHamiltonian{true}, x, p) = (sys.H_p(buf, x, p); buf)
@inline _eval_Hp!(buf, sys::FreidlinWentzellHamiltonian{false}, x, p) = copyto!(buf, sys.H_p(x, p))

raw"""
In-place evaluation of ``\partial_x H``. See [`_eval_Hp!`](@ref).
"""
@inline _eval_Hx!(buf, sys::FreidlinWentzellHamiltonian{true}, x, p) = (sys.H_x(buf, x, p); buf)
@inline _eval_Hx!(buf, sys::FreidlinWentzellHamiltonian{false}, x, p) = copyto!(buf, sys.H_x(x, p))

function FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    ds isa CoupledSDEs && proper_FW_system(ds)
    a = _trace_normalized_a(ds)
    f = dynamic_rule(ds)
    ps = current_parameters(ds)
    iip = Val(SciMLBase.isinplace(ds))
    H_p = _make_H_p(a, f, ps, iip)
    H_x = _make_H_x(a, ds, f, ps, iip)
    x_ref = collect(current_state(ds))
    return FreidlinWentzellHamiltonian{
        true, D, typeof(H_x), typeof(H_p), typeof(a), typeof(x_ref),
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

_make_H_p(a, f, ps, ::Val{false}) = _AutoHamiltonianPOOP(a, f, ps)
_make_H_p(a, f, ps, ::Val{true}) = _AutoHamiltonianPIIP(a, f, ps)

function _make_H_x(a, _ds, f, ps, ::Val{false})
    jac_fn = let f = f, ps = ps
        x -> f(x, ps, 0.0)
    end
    return _AutoHamiltonianXOOP(a, jac_fn)
end

_make_H_x(a, ds, _f, ps, ::Val{true}) = _AutoHamiltonianXJac(a, jacobian(ds), ps)

_add_da_term!(::Val{true}, _, _, _, _, _) = nothing

function _add_da_term!(::Val{false}, out, a, x, p, x_probe)
    h = _fd_step(eltype(x))
    inv_2h = inv(2 * h)
    copyto!(x_probe, x)
    @inbounds for l in eachindex(out)
        x_probe[l] = x[l] + h
        a_xp = a(x_probe)
        x_probe[l] = x[l] - h
        a_xm = a(x_probe)
        x_probe[l] = x[l]
        contraction = zero(eltype(out))
        for j in eachindex(p), k in eachindex(p)
            contraction += p[j] * (a_xp[j, k] - a_xm[j, k]) * p[k]
        end
        out[l] += 0.5 * contraction * inv_2h
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
