using CriticalTransitions, StaticArrays
CT = CriticalTransitions

function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x-x^3 -10*x*y^2
    dy = -(1+x^2)*y
    SA[dx, dy]
end
σ = 0.25

sys = StochSystem(meier_stein, [], zeros(2))

tspan = (0.0, 100.0)
prob = ODEProblem{false}(meier_stein, SA[0.0, 0.0], tspan, ())
kpo = CoupledODEs(prob) # DynamicalSystems method

using Test

# Test StochSystem constructors
@testset "StochSystem constructors" begin
    f(x,p,t) = x^2
    sys1 = StochSystem(f, [1], [0])
    @test sys1.f(2, sys1.pf, 0) == 4
    @test sys1.u == [0]
    @test sys1.σ == 0.0
    @test sys1.g == idfunc
    @test sys1.pg == nothing
    @test sys1.Σ == I(1)
    @test sys1.process == "WhiteGauss"

    sys2 = StochSystem(f, [1], [0], 0.5)
    @test sys2.f(2, sys2.pf, 0) == 4
    @test sys2.u == [0]
    @test sys2.σ == 0.5
    @test sys2.g == idfunc
    @test sys2.pg == nothing
    @test sys2.Σ == I(1)
    @test sys2.process == "WhiteGauss"

    Σ = [true false; false true]
    sys3 = StochSystem(f, [1], [0], 0.5, Σ)
    @test sys3.f(2, sys3.pf, 0) == 4
    @test sys3.u == [0]
    @test sys3.σ == 0.5
    @test sys3.g == idfunc
    @test sys3.pg == nothing
    @test sys3.Σ == Σ
    @test sys3.process == "WhiteGauss"
end

# Test drift function
@testset "drift" begin
    f(x,p,t) = x.^2
    sys = StochSystem(f, [1], [0])
    @test drift(sys, [2]) == [4]
end

# Test σg function
@testset "σg" begin
    f(x,p,t) = x^2
    sys = StochSystem(f, [1], [0], 0.5)
    @test CT.σg(sys)(2, sys.pf, 0) == 0.5 .* CT.idfunc(2, sys.pf, 0)
end

# Test p function
@testset "p" begin
    f(x,p,t) = x^2
    sys = StochSystem(f, [1], [0], 0.5)
    @test CT.p(sys) == [[1], nothing]
end

# Test CoupledODEs function
@testset "CoupledODEs" begin
    f(x,p,t) = x.^2
    sys = StochSystem(f, [1.0], [0.0], 0.5)
    cds = CoupledODEs(sys)
    @test cds.integ.f([2], sys.pf, 0) == [4]
    @test cds.integ.u == [0]
    @test cds.p0 == sys.pf
end


# Test StochSystem conversion from CoupledODEs
@testset "StochSystem from CoupledODEs" begin
    f(x,p,t) = x.^2
    cds = CoupledODEs(f, [0.0], [0.0])
    sys = StochSystem(cds, 0.5, idfunc, nothing, I(1), "WhiteGauss")
    @test sys.f([2], sys.pf, 0) == [4]
    @test sys.u == [0]
    @test sys.σ == 0.5
    @test sys.g == idfunc
    @test sys.pg == nothing
    @test sys.Σ == I(1)
    @test sys.process == "WhiteGauss"
end
