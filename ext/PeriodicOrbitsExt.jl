module PeriodicOrbitsExt

using CriticalTransitions: CriticalTransitions, LimitCycleFrame, CoupledSDEs
using PeriodicOrbits: PeriodicOrbit

function CriticalTransitions.LimitCycleFrame(po::PeriodicOrbit, sys::CoupledSDEs; kw...)
    return LimitCycleFrame(po.points, po.T, sys; kw...)
end

end # module
