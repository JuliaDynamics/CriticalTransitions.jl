function initial_curve()
    println("in initial_curve()")

    # Read initial curve points from file
    if !isfile(ficurve_name)
        error("Please provide a data file for the initial curve")
    end

    # Initialize curve from file
    points = Vector{Tuple{Float64,Float64}}()
    open(ficurve_name, "r") do f
        for line in eachline(f)
            x, y = parse.(Float64, split(line))
            push!(points, (x, y))
        end
    end

    n = length(points)
    if n > NCURVEMAX
        error("Error: j = $n while NCURVEMAX = $NCURVEMAX")
    end

    # Create linked list of curve points
    global NCURVE = n
    global icurve = Vector{MyCurve}(undef, n)

    for i in 1:n
        x, y = points[i]
        next_idx = i == n ? 1 : i + 1
        prev_idx = i == 1 ? n : i - 1
        icurve[i] = MyCurve(x, y, nothing, nothing)
    end

    # Link the nodes
    for i in 1:n
        icurve[i] = MyCurve(
            icurve[i].x,
            icurve[i].y,
            i == n ? icurve[1] : icurve[i+1],
            i == 1 ? icurve[n] : icurve[i-1]
        )
    end

    println("start")
    Nac = 0
    iac = zeros(Int, NX*100)

    # Process each curve segment
    for n in 1:NCURVE
        xc0, yc0 = icurve[n].x, icurve[n].y
        next_point = icurve[n].next
        xc1, yc1 = next_point.x, next_point.y

        i0 = min(floor(Int, (xc0 - XMIN)/hx), floor(Int, (xc1 - XMIN)/hx))
        j0 = min(floor(Int, (yc0 - YMIN)/hy), floor(Int, (yc1 - YMIN)/hy))
        i1 = max(ceil(Int, (xc0 - XMIN)/hx), ceil(Int, (xc1 - XMIN)/hx))
        j1 = max(ceil(Int, (yc0 - YMIN)/hy), ceil(Int, (yc1 - YMIN)/hy))

        for i in i0:i1, j in j0:j1
            ind = i + j*NX
            if ms[ind] == 0
                g[ind] = init(getpoint(ind))
                ms[ind] = 1
                solinfo[ind] = SolInfo('0', 0, 0, 0.0)
                addtree(ind)
                Nac += 1
                iac[Nac] = ind
            end
        end
    end
end

function init(x::MyVector)::Float64
    if chfield == 'l'
        x = vec_difference(x, x_ipoint)
        a, b = -2.0, -alin
        c, d = 2.0*alin, -1.0
        aux1 = c - b
        aux2 = a + d
        aux = aux1*aux1 + aux2*aux2
        aux1 *= aux2/aux
        aux2 *= aux2/aux
        A = -(a*aux2 + c*aux1)
        B = -(b*aux2 + d*aux1)
        C = -(d*aux2 - b*aux1)
        return A*x.x*x.x + 2.0*B*x.x*x.y + C*x.y*x.y

    elseif chfield == 'c'
        # Find closest point on curve
        step = max(1, NCURVE÷12)
        nmin, dmin = 0, INFTY

        # Course search
        for n in 1:step:NCURVE
            dtemp = length(icurve[n].x - x.x, icurve[n].y - x.y)
            if dtemp < dmin
                dmin = dtemp
                nmin = n
            end
        end

        # Fine search forward
        imin = nmin
        step *= 2
        curr = icurve[nmin]
        for _ in 1:step
            curr = curr.next
            dtemp = length(curr.x - x.x, curr.y - x.y)
            if dtemp < dmin
                dmin = dtemp
                imin = (nmin + n) % NCURVE
            end
        end

        # Fine search backward
        curr = icurve[nmin]
        for _ in 1:step
            curr = curr.prev
            dtemp = length(curr.x - x.x, curr.y - x.y)
            if dtemp < dmin
                dmin = dtemp
                imin = (nmin + NCURVE - n) % NCURVE
            end
        end

        # Calculate result
        x0 = MyVector(icurve[imin].x, icurve[imin].y)
        x1 = MyVector(icurve[imin].next.x, icurve[imin].next.y)
        v0 = vec_difference(x0, x)
        v1 = vec_difference(x1, x)

        if dot_product(v0, vec_difference(v0, v1)) < 0.0
            x1 = MyVector(icurve[imin].prev.x, icurve[imin].prev.y)
            v1 = vec_difference(x1, x)
        end

        v = vec_difference(x1, x0)
        c = length_vec(v)
        dmin = abs(v0.x*v1.y - v0.y*v1.x)/c
        a = sqrt(v0.x*v0.x + v0.y*v0.y - dmin*dmin)
        aux = a/c
        v = MyVector(v.x * aux, v.y * aux)
        x0 = vec_sum(x0, v)
        x1 = vec_lin_comb(x, x0, 0.5, 0.5)
        bvec0 = myfield(x0)
        bvec1 = myfield(x1)
        bvec = myfield(x)
        lb0 = length_vec(bvec0)
        aux = dot_product(bvec0, bvec1)/(lb0*lb0)
        bvec0.x *= aux
        bvec0.y *= aux
        aux = dot_product(bvec0, bvec)/(lb0*lb0)
        bvec.x *= aux
        bvec.y *= aux
        temp = dmin*(4.0*length_vec(vec_difference(bvec1, bvec0)) + length_vec(vec_difference(bvec, bvec0)))/3.0
        return temp
    else
        return 0.0
    end
end
