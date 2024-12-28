function string_method(
    sys::SgmamSystem,
    x_initial;
    ϵ::Float64=1e-1,
    iterations::Int64=1000,
    show_progress::Bool=false,
)
    H_p, H_x = sys.H_p, sys.H_x
    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, alpha = init_allocation_string(x_initial, Nt)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        x += ϵ*H_p(x,0*x)
        # reset initial and final points to allow for string computation
        # between points that are not stable fixed points
        x[:,1] = x_initial[:,1]; x[:,end] = x_initial[:,end];

        # reparametrize to arclength
        interpolate_path!(x, alpha, s)

        next!(progress; showvalues=[("iterations", i),])
    end
    return StateSpaceSet(x)
end
function init_allocation_string(x_initial, Nt)
    # preallocate
    x = deepcopy(x_initial)
    alpha = zeros(Nt)
    return x, alpha
end

function string_method(
    sys::CoupledSDEs,
    x_initial;
    ϵ::Float64=1e-1,
    iterations::Int64=1000,
    show_progress::Bool=false,
)
    b(x) = drift(sys, x)
    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, alpha = init_allocation_string(x_initial, Nt)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        # do not touch the initial and final points to allow for string computation
        # between points that are not stable fixed points
        for j in 2:(size(x)[2]-1)
            x[:,j] += ϵ*b(x[:,j])
        end

        # reparametrize to arclength
        interpolate_path!(x, alpha, s)

        next!(progress; showvalues=[("iterations", i),])
    end
    return StateSpaceSet(x)
end
