"""
$(TYPEDEF)

A struct representing a Langevin dynamical system with damping rate `gamma`` and temperature `beta``.

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
    "Function computing the kinetic energy of the system."
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
