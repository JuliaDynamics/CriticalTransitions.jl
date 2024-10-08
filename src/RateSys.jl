#  |                    |                        |                  |
# t_i    autonomous    t_ni   non-autonomous   t_nf   autonomous   t_f
#  |                    |                        |                  |

using DynamicalSystems
using Plots

# we write

function RateSyst(tni,tnf,f,λ,p_λ,initvals)
    func(u,p,t) = combined_system(u,t,tni,tnf,f,λ,p_λ);
    return ContinuousDynamicalSystem(func, initvals, Float64[], t0=tni)
end

function combined_system(u,t,tni,tnf,f,λ,p_λ)
    lambda = t < tni ? λ(u,p_λ,tni) : tni <= t <= tnf ? λ(u,p_λ,t) : λ(u,p_λ,tnf)
    return f(u,lambda,t)
end;


###############
# user writes

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
sys = RateSyst(tni,tnf,f,λ,p_λ,initvals)

traj=trajectory(sys,40.,t0=-20.)

trajic=[traj[1][i][1] for i in 1:length(traj[1])]
plot(collect(traj[2]),trajic)


# further ideas
# function futureSyst(RateSyst)
# Canard Trajectories

# \dot(x) = f(x(t),λ(t))



