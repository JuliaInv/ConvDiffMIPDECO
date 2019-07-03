using ConvDiffMIPDECO
using Test
using jInv.Mesh
using MAT
using PyPlot


domain = [0. 3. 0 1. 0 2.]
n      = 3*[7 9 12]-1
M      = getRegularMesh(domain,n)
Mass, Mass_const, SM = getFEMMatrices3D(M)

e = ones(prod(M.n+1))
@test abs(prod((domain[2:2:end]-domain[1:2:end])) - dot(e,Mass*e))/dot(e,Mass*e) < 1e-2

f = getFEMsource3D(M)
v = Mass_const*f
u = SM\v