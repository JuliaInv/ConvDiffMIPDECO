using ConvDiffMIPDECO
using Test
using jInv.Mesh


domain = [0. 3. 0 1.]
n      = [60 20].-1
M      = getRegularMesh(domain,n)
@time Mass, Mass_const, SM = getFEMMatrices2D(M)


f = getFEMsource2D(M)
v = Mass_const*f
u = SM\v