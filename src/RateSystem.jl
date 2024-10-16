#  |                    |                        |                  |
# t_i    autonomous    t_ni   non-autonomous   t_nf   autonomous   t_f
#  |                    |                        |                  |

# we write

function RateSystem(tni,tnf,f,λ,p_λ,initvals)
    func(u,p,t) = combined_system(u,t,tni,tnf,f,λ,p_λ);
    return ContinuousTimeDynamicalSystem(func, initvals, Float64[], t0=tni)
end

function combined_system(u,t,tni,tnf,f,λ,p_λ)
    lambda = t < tni ? λ(u,p_λ,tni) : tni <= t <= tnf ? λ(u,p_λ,t) : λ(u,p_λ,tnf)
    return f(u,lambda,t)
end;


###############
# user writes

# ...moved to test/RateSystem.jl (Reyk)


# further ideas
# function futureSyst(RateSyst)
# Canard Trajectories

# \dot(x) = f(x(t),λ(t))



