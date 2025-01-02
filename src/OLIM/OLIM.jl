module OLIM

const INFTY = Inf
const TOL = 1.0e-12
const NX = 1024
const NY = 1024
const K = 22
const NCURVEMAX = 399

const nx1 = NX - 1
const ny1 = NY - 1
const nxy = NX * NY

const KK = K * K

NCURVE = 0
count = 0
h = 0.0
hx = 0.0
hy = 0.0

ms = fill(0, NX * NY)
g = fill(INFTY, NX * NY)
pos = fill(0, NX * NY)
tree = fill(0, NX * NY)

const neii = [1, NX + 1, NX, NX - 1, -1, -NX - 1, -NX, -NX + 1]


Uexact = fill(0.0, NX * NY)
solinfo = [SolInfo('n', 0, 0, 0.0) for _ in 1:NX*NY]
NXY = NX * NY


NMESH = 0
x_ipoint = [0.0, 0.0]
icurve = Vector{MyCurve}(undef, NCURVEMAX)


# Variables for the potential
# // define the computational domain
XMIN = 0.0
XMAX = 0.0
YMIN = 0.0
YMAX = 0.0





include("types.jl")
include("point.jl")
include("field.jl")


end # module OLIM
