# CT.jl Vision

## Main types

* DynamicalSystem
    * StochSystem(f, pf, u, g, pg, sigma, Sigma, process)
    * RateSystem(f, pf, u, r(pf, t))


* Methods
    * simulate
    * relax
    * to_ds

lyapunovspectrum(CoupledODEs(sys))

    r(pf, pr, t)
    pr = T_tr, T_ramp, ...



    function f(u,p,t)
        a, b = p
        x, y = u

        dx = ...
        dy = ...

        return [dx, dy]


    function r(p, t)
        a, b = p

        da = ...
        db = ...


    g(u,p,t) = f(u, r(p,t), t)

    relax(rsys, t)
        ODEProblem(rsys.f, rsys.r, rsys.u)

    relax(ssys)
        ODEproble