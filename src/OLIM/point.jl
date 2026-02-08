

function getpoint(ind)
    if ind < NXY
        x = hx * (mod(ind, NX)) + XMIN
        y = hy * (div(ind, NX)) + YMIN
        return (x=x, y=y)  # Returns a NamedTuple as a vector-like structure
    end
end


function ipoint()
    # Array for surface indices offset
    isur = [0, 1, NX+1, NX]

    # Calculate initial indices
    i = floor(Int, (x_ipoint.x - XMIN)/hx)
    j = floor(Int, (x_ipoint.y - YMIN)/hy)
    ind0 = i + j*NX

    # Process the four surrounding points
    for m in 1:4
        ind = ind0 + isur[m]
        x = getpoint(ind)
        gtemp = init(x)
        g[ind] = gtemp
        ms[ind] = 1
        addtree(ind)
        solinfo[ind].type = 0
    end
end
