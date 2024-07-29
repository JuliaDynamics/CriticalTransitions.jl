using Symbolics,OptControl

@variables t u[1:2] x[1:2]
f = [0 1; 0 0] * x + [0, 1] * u
L = 0.5 * u^2
t0 = [1.0, 1.0]
tf = [0.0, 0.0]
tspan = (0.0, 2.0)
N = 100
sol = generateJuMPcodes(L, f, x, u, tspan, t0, tf; N=N, writeFilePath="examples/JuMPscript.jl")
