## the following function integrates the nonautonomous system created by the ContinuousTimeDynamicalSystem ds and RateProtocol rp, with rate-parameter r 
## the initial condition is e_start at time t_start
## the system is integrated for T = t_end-t_start time units and the trajectory position at time T is then compared with the location of e_end 
## if the trajectory position at time T is within a neighbourhood of radius rad of e_end, the function returns true, and false otherwise  

function end_point_tracking(
    ds::ContinuousTimeDynamicalSystem, rp::RateProtocol, r::Float64, e_start::Vector, t_start::Float64, e_end::Vector, t_end::Float64, rad::Float64
)

    # note that should ensure that Δt << dλ(rt)/dt
    # more generally, the suitability of the numerical integrator diffeq has to be ensured

    T = t_end - t_start
    rp.r = r # set the rate-parameter to the minimum value 
    rate_sys_min_rate = apply_ramping(ds, rp, t_start) # create the nonautonomous system
    traj_min_rate = trajectory(rate_sys_min_rate, T, e_start) # simulate the nonautonomous system

    return abs(traj_min_rate[1][end,:]-e_end) < rad 

end

# the following function finds via brute force an interval [a,b] satisfying b-a ≤ tol satisfying f(a)*f(b) == false for a function f returning boolean values 
# i.e. a small interval [a,b] for which the value of the function f evaluated at the points a and b returns both true and false  

function binary_bisection(f::Function, a::Float64, b::Float64;
    tol::Float64=1e-12, maxiter::Int64=10000)
    a < b || error("Please enter a, b such that a < b.") 
    fa = f(a)
    fa*f(b) == false || error("The function does not satisfy the initial requirement that f(a)*f(b) == false")
    ii = 0
    local c
    while b-a > tol
        ii += 1
        ii ≤ maxiter || error("The maximum number of iterations was exceeded.")
        c = (a+b)/2
        fc = f(c)
        if fa*fc == true # reduce to the shrunken interval [c,b], since we must have f(c)*f(b) == false
            a = c  
            fa = fc 
        else # reduce to the shrunken interval [a,c], since we have f(a)*f(c) == false
            b = c  
        end
    end
    return [a,b]
end;

## the following function aims to find a critical value for the rate-parameter between successful and non-successful "end-point tracking" 

## the end-point tracking is currently defined with respect to a stable equilibria of the frozen-in system attained at time t = t_start...
##...and a stable equilibria of the frozen-in system attained at time t = t_end 

## there is currently no check that e_start and e_end are connected via a continuous branch of moving sinks attained across [t_start,t_end],...
##... which would be a desirable check to perform although requires more work to implement (for a first draft I skipped this)
## before doing this it probably makes sense to modify the moving_sinks function so as to return continuous branches of grouped equilibria...
##...rather than just a vector of vectors containing all the (non-grouped) equilibria found for a range of frozen-in times  

## there are check for ensuring that e_start and e_end are indeed stable equilibria of the frozen-in systems attained at times t_start and t_end respectively
## there are also checks that the system succesfully end-point tracks for r = rmin and fails to end-point track for r = rmax (or vice versa, for exotic cases)

## the function ultimately passes the end_point_tracking function to the binary_bisection function to find a critical rate for end-point tracking
 

function critical_end_point_tracking_rate(
    ds::ContinuousTimeDynamicalSystem, rp::RateProtocol, e_start::Vector, t_start::Float64, e_end::Vector, t_end::Float64;  
    rmin::Float64 = 1.e-2, 
    rmax::Float64 = 1.e2,
    rad::Float64 = 1.e-1,
    tol::Float64 = 1.e-3,
    maxiter::Int64 = Int64(1.e4)
)

    rmin < rmax || error("Please enter rmin, rmax such that rmin < rmax.") 

    ##-----check no.1-----## 

    ## that e_start is a stable equilibrium of the frozen-in system attained at t = t_start

    rate_sys_start = apply_ramping(ds, rp, t_start)
    box_start = [interval(x-0.1,x+0.1) for x ∈ e_start]
    fp_start, ~, stab_start = fixedpoints(rate_sys_start, box_start)
    if any(e -> isapprox(e,e_start;atol=1e-4),fp_start)
        idx = findmin([abs(e-e_start) for e ∈ fp_start])[2] # the index in fp_start of the fixed point which is approximately equal to e_end
        stab_start[idx] || @warn("The vector e_start = $(e_start) does not describe a stable equilibrium of the frozen-in system attained at t = t_start = $(t_start).")
    else
        @warn("The vector e_start = $(e_start) does not describe a stable equilibrium of the frozen-in system attained at t = t_start = $(t_start).")
    end

    ##-----check no.2-----##
    
    ## that e_end is a stable equilibrium of the frozen-in system attained at t = t_end

    rate_sys_end = apply_ramping(ds, rp, t_end)
    box_end = [interval(x-0.1,x+0.1) for x ∈ e_end]
    fp_end, ~, stab_end = fixedpoints(rate_sys_end, box_end)
    if any(e -> isapprox(e,e_end;atol=1e-4),fp_end)
        idx = findmin([abs(e-e_end) for e ∈ fp_end])[2] # the index in fp_end of the fixed point which is approximately equal to e_end
        stab_end[idx] || @warn("The vector e_end = $(e_end) does not describe a stable equilibrium of the frozen-in system attained at t = t_end = $(t_end).")
    else
        @warn("The vector e_end = $(e_end) does not describe a stable equilibrium of the frozen-in system attained at t = t_end = $(t_end).")
    end

    ##-----check no.3-----##
    
    ## that
    ## the system successfully end-point tracks for r = rmin and fails to end-point track for r = rmax
    ## or
    ## the system fails to end-point track for r = rmin and successfully end point tracks for r = rmax 

    min_rate_track = end_point_tracking(ds,rp,rmin,e_start,t_start,e_end,t_end,rad) # true if the system end-point tracks for that rate

    max_rate_track = end_point_tracking(ds,rp,rmax,e_start,t_start,e_end,t_end,rad) # true if the system end-point tracks for that rate

    if min_rate_track && max_rate_track 
        error("The system successfully end-point tracks for both r = rmin = $(rmin) and r = rmax = $(rmax).")
    elseif ~min_rate_track && ~max_rate_track 
        error("The system fails to end-point track for both r = rmin = $(rmin) and r = rmax = $(rmax).")
    elseif ~min_rate_track && max_rate_track 
        @warn("The system fails to end-point track for r = rmin = $(rmin) but succesfully end-point tracks for r = rmax = $(rmax).")
    end

    ## all checks complete ##

    ## performing the bisection to find a critical rate within the given tolerance ##

    func(r) = end_point_tracking(ds,rp,r,e_start,t_start,e_end,t_end,rad)
    rcrit_interval = binary_bisection(func,rmin,rmax;tol,maxiter)

    return (rcrit_interval[1]+rcrit_interval[2])/2

end

# beginning notes

# you are trying to locate the critical rate r for which the system end-point tracks the moving branch of sinks 

# the first case is where after some travel time you are within a defined neighbourhood of a sink belonging to the frozen-in system at that time 
# for this you should provide the tuple of tuples ((e₋,t₋),(e₊,t₊))
# then you start at the position e⁻ at the time t⁻ and simulate t⁺-t⁻ time-units and compute the distance to e⁺
# we assume that there exists some pairing (r₁,r₂) whereby the system lies within the neighbourhood of e⁺ under the rate r₁ and outside of the neighbourhood of e⁺ under the rate r₂ 
