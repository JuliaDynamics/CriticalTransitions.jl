



# Main function converted to Julia style
function main()
    param()
    if chfield in ('l', 'q', 'r')
        ipoint()
    elseif chfield in ('c', 'b')
        initial_curve()
    else
        error("Invalid chfield: $chfield")
    end

    t_start = time()
    olim()
    cpu_time = time() - t_start
    println("cputime of olim() = $cpu_time")

    # Write results to files
    open(f_qpot_name, "w") do fg
        open(f_solinfo_name, "w") do fs
            k = 0
            errmax = 0.0
            erms = 0.0

            for j in 1:NY, i in 1:NX
                ind = i + NX * (j - 1)
                # ...existing output logic...
            end
        end
    end

    println("NX = $NX, NY = $NY")
    if chfield in ('l', 'c')
        println("errmax = $errmax, erms = $(sqrt(erms/k))")
    end
end

main()
