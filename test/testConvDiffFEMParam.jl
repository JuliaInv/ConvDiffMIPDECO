using Test
using ConvDiffMIPDECO
using jInv.Mesh
using jInv.Utils
using jInv.ForwardShare
using jInv.LinearSolvers


println("\t\t-- test ConvDiffFEMParam --")
domain = [-1 1. -1 1]
n      = [32 32]
M      = getRegularMesh(domain,n)
Ainv   = getJuliaSolver()
Ainv   = ConvDiffMIPDECO.getBICGSTB(PC=:ssor)
Ainv   = getMUMPSsolver()
pFor   = getConvDiffFEMParam(M,sig=.0001,Ainv=Ainv)
m0     = zeros(tuple(n...))
m0[10:20,15:25] =1.;
isOK, his = checkDerivative(vec(m0),pFor,out=true)
@test isOK

v = randn(M.nc)
w = randn(prod(M.n.+1))
t1 = dot(w, getSensMatVec(v,vec(m0),pFor))
t2 = dot(v, getSensTMatVec(w,vec(m0),pFor))
@test abs(t1-t2)/abs(t1) < 1e-8

domain = [-1 1. 0 1 -1 1]
n      = [6 8 10]
M      = getRegularMesh(domain,n)
clear!(Ainv)
pFor   = getConvDiffFEMParam(M,sig=.5,Ainv=Ainv)

m0              = zeros(tuple(n...))
m0[3:4,4:6,4:7] =1.;
isOK, his = checkDerivative(vec(m0),pFor,out=true)
@test isOK

v = randn(M.nc)
w = randn(prod(M.n.+1))
t1 = dot(w, getSensMatVec(v,vec(m0),pFor))
t2 = dot(v, getSensTMatVec(w,vec(m0),pFor))
@test abs(t1-t2)/abs(t1) < 1e-8

