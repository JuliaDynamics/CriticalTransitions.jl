using Plots

function f(u::SVector{1, Float64},p::SVector{1, Float64},t::Float64)
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;

function λ(p::Vector{Float64},t::Float64)
    r,λ_max = p
    lambda = (λ_max/2)*(tanh(λ_max*r*t/2) +1)
    return SVector{1}(lambda)
end;

tni=-10.; tnf=10.; p_λ = [1.0,3.0]; initvals = [-4.];
sys = RateSystem(tni,tnf,f,λ,p_λ,initvals)

traj=trajectory(sys,40.,t0=-20.)

trajic=[traj[1][i][1] for i in 1:length(traj[1])]
plot(collect(traj[2]),trajic)