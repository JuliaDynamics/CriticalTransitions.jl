# Test for idfunc
@testset "idfunc" begin
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    expected = [1, 1, 1]
    @test idfunc(u, p, t) == expected
end

# Test for idfunc!
@testset "idfunc!" begin
    du = zeros(3)
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    expected = [1, 1, 1]
    idfunc!(du, u, p, t)
    @test du == expected
end

# Test for is_iip
@testset "is_iip" begin
    @test is_iip(idfunc) == false
    @test is_iip(idfunc!) == true
end

# Test for intervals_to_box
@testset "intervals_to_box" begin
    using CriticalTransitions.IntervalArithmetic
    bmin = [-2, -1, 0]
    bmax = [2, 1, 1]
    expected = (Interval([-2, 2]...), Interval([-1, 1]...), Interval([0, 1]...))
    @test intervals_to_box(bmin, bmax).v.data == expected
end

# Test for additive_idx!
@testset "additive_idx!" begin
    du = zeros(3)
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    idx = 2
    expected = [0, 1, 0]
    additive_idx!(du, u, p, t, idx)
    @test du == expected
end

# Test for additive_idx
@testset "additive_idx" begin
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    idx = 2
    expected = SVector{3}([0, 1, 0])
    @test additive_idx(u, p, t, idx) == expected
end

# Test for multiplicative_idx!
@testset "multiplicative_idx!" begin
    du = zeros(3)
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    idx = 2
    expected = [0, 2, 0]
    multiplicative_idx!(du, u, p, t, idx)
    @test du == expected
end

# Test for multiplicative_idx
@testset "multiplicative_idx" begin
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    idx = 2
    expected = SVector{3}([0, 2, 0])
    @test multiplicative_idx(u, p, t, idx) == expected
end

# Test for multiplicative_first!
@testset "multiplicative_first!" begin
    du = zeros(3)
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    expected = [1, 0, 0]
    CT.multiplicative_first!(du, u, p, t)
    @test du == expected
end

# Test for multiplicative_first
@testset "multiplicative_first" begin
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    expected = SVector{3}([1, 0, 0])
    @test CT.multiplicative_first(u, p, t) == expected
end

# Test for additive_first!
@testset "additive_first!" begin
    du = zeros(3)
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    expected = [1, 0, 0]
    CT.additive_first!(du, u, p, t)
    @test du == expected
end

# Test for additive_first
@testset "additive_first" begin
    u = [1, 2, 3]
    p = [0.1, 0.2, 0.3]
    t = 0.5
    expected = SVector{3}([1, 0, 0])
    @test CT.additive_first(u, p, t) == expected
end

# Test for anorm
@testset "anorm" begin
    vector = [1, 2, 3]
    A = [1 0 0; 0 2 0; 0 0 3]
    expected = sqrt(1^2 + 2^2 * 2 + 3^2 * 3)
    @test anorm(vector, A) == expected
end

# Test for subnorm
@testset "subnorm" begin
    vector = [3, 7, 4]
    directions = [1, 3]
    expected = sqrt(3^2 + 4^2)
    @test subnorm(vector; directions=directions) == expected
end

# Test for central
@testset "central" begin
    f = [1, 2, 3, 4, 5]
    idx = 3
    dz = 1
    expected = (f[idx + 1] - f[idx - 1]) / (2 * dz)
    @test CT.central(f, idx, dz) == expected
end

# Test for central2
@testset "central2" begin
    f = [1, 2, 3, 4, 5]
    idx = 3
    dz = 1
    expected = (f[idx + 1] - 2 * f[idx] + f[idx - 1]) / dz^2
    @test CT.central2(f, idx, dz) == expected
end

# Test for smoothabs
@testset "smoothabs" begin
    x = 2
    xi = 1000
    expected = x * tanh(x * xi)
    @test CT.smoothabs(x, xi) == expected
end
