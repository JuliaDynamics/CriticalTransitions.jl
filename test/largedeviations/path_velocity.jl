function reference_path_velocity(path, time; order = 4)
    v = similar(path)
    N = size(path, 2)
    if order == 2
        @inbounds @views begin
            inv_h1 = 1 / (time[2] - time[1])
            inv_hN = 1 / (time[end] - time[end - 1])
            @. v[:, 1] = (path[:, 2] - path[:, 1]) * inv_h1
            @. v[:, end] = (path[:, end] - path[:, end - 1]) * inv_hN
            for i in 2:(N - 1)
                inv_hi = 1 / (time[i + 1] - time[i - 1])
                @. v[:, i] = (path[:, i + 1] - path[:, i - 1]) * inv_hi
            end
        end
    elseif order == 4
        @inbounds @views begin
            inv_h1 = 1 / (time[2] - time[1])
            inv_hN = 1 / (time[end] - time[end - 1])
            inv_h2 = 1 / (time[3] - time[1])
            inv_hM = 1 / (time[end] - time[end - 2])
            @. v[:, 1] = (path[:, 2] - path[:, 1]) * inv_h1
            @. v[:, end] = (path[:, end] - path[:, end - 1]) * inv_hN
            @. v[:, 2] = (path[:, 3] - path[:, 1]) * inv_h2
            @. v[:, end - 1] = (path[:, end] - path[:, end - 2]) * inv_hM
            for i in 3:(N - 2)
                inv6 = 1 / (6 * (time[i + 1] - time[i - 1]))
                @. v[:, i] = (
                    -path[:, i + 2] + 8 * path[:, i + 1] -
                        8 * path[:, i - 1] + path[:, i - 2]
                ) * inv6
            end
        end
    end
    return v
end

@testset "path velocity scalar stencils" begin
    time = [0.0, 0.07, 0.18, 0.31, 0.47, 0.66, 0.82, 1.0]
    path = [
        sin.(time)';
        cos.(2 .* time)';
        (time .^ 3 .- 0.4 .* time)';
    ]

    for order in (2, 4)
        expected = reference_path_velocity(path, time; order)
        actual = fill(NaN, size(path))
        @test CT.path_velocity!(actual, path, time; order) === actual
        @test actual ≈ expected rtol = 0 atol = 8eps(Float64)
        @test CT.path_velocity(path, time; order) ≈ expected rtol = 0 atol = 8eps(Float64)
    end

    untouched = fill(3.0, size(path))
    CT.path_velocity!(untouched, path, time; order = 3)
    @test all(==(3.0), untouched)
end
