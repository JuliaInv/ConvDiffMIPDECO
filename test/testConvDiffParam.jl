using Test
using ConvDiff
using jInv.Mesh
using jInv.Utils
using jInv.ForwardShare

println("\t\t-- test ConvDiffParam --")
domain = [-1 1. -1 1]
n      = [32 32]
M      = getRegularMesh(domain,n)
bc     = (:dir,:neu,:dir,:neu)
h      = M.h
v      = X -> [ones(size(X,1)) zeros(size(X,1))]
gd     = X -> X[:,1].*X[:,2]
gn     = X -> ones(size(X,1))
pFor   = getConvDiffParam(M,v,sig=.5,gd=gd,gn=gn,bc=bc)

m0     = zeros(tuple(n...))
m0[10:20,15:25] =1.;
isOK, his = checkDerivative(vec(m0),pFor,out=true)
@test isOK
v = randn(M.nc)
w = randn(prod(M.nc))
t1 = dot(w, getSensMatVec(v,vec(m0),pFor))
t2 = dot(v, getSensTMatVec(w,vec(m0),pFor))
@test abs(t1-t2)/abs(t1) < 1e-8

domain = [-1 1. 0 1 -1 1]
n      = [6 8 10]
M      = getRegularMesh(domain,n)
bc     = (:dir,:neu,:dir,:neu,:dir,:neu)
h      = M.h
v      = X -> [ones(size(X,1)) zeros(size(X,1)) zeros(size(X,1))]
gd     = X -> X[:,1].*X[:,2]
gn     = X -> ones(size(X,1))
pFor   = getConvDiffParam(M,v,sig=.5,bc=bc)

m0              = zeros(tuple(n...))
m0[3:4,4:6,4:7] =1.;
isOK, his = checkDerivative(vec(m0),pFor,out=true)
@test isOK
v = randn(M.nc)
w = randn(prod(M.nc))
t1 = dot(w, getSensMatVec(v,vec(m0),pFor))
t2 = dot(v, getSensTMatVec(w,vec(m0),pFor))
@test abs(t1-t2)/abs(t1) < 1e-8
