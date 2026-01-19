"""
$(TYPEDEF)

A stochastic dynamical system obeying Langevin dynamics of the form.
```math
\\dot x = p \\,, \\\\
\\dot p = - \\gamma p - \\nabla U(x) + \\sqrt{2\\gamma\\beta^{-1}} \\dot W_t \\,,
```
with damping coefficient ``\\gamma`` and inverse temperature ``\\beta``. The Hamiltonian
``H = U + K`` is given by the potential energy ``U`` and kinetic energy ``K = p^2/2``.

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

"""
struct LangevinSystem{H,D,KE,T}
    "Function that computes the total energy (kinetic + potential) of the system."
    Hamiltonian::H
    "Function representing the divergence-free part of the drift."
    driftfree::D
    "Function giving the kinetic energy of the system."
    kinetic::KE
    "Damping coefficient that determines the strength of friction."
    gamma::T
    "Inverse temperature parameter (β = 1/kT) that sets the noise intensity."
    beta::T
end

function SciMLBase.remake(
    sys::LangevinSystem;
    Hamiltonian=sys.Hamiltonian,
    driftfree=sys.driftfree,
    kinetic=sys.kinetic,
    gamma=sys.gamma,
    beta=sys.beta,
)
    return LangevinSystem(Hamiltonian, driftfree, kinetic, gamma, beta)
end
