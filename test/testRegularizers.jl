using ConvDiff
using Base.Test
using jInv.Mesh
using jInv.Utils

domain = [1.2 4.5 1.2 2.4]
n      = [4 7]
M      = getRegularMesh(domain,n)

mc     = rand(M.nc)

isOK, = checkDerivative(m->mipReg(m,0*m,M)[1], m->(mipReg(m,0*m,M)[2])',mc)
@test isOK

