using jInv.LinearSolvers
using Test
using ConvDiff
using jInv.ForwardShare
using jInv.Mesh
using LinearAlgebra
using SparseArrays

domain = [0 1. 0 1]

N = 2*[16;32;64;128]
v = randn(2)
v = v/norm(v)
V = X -> [v[1]*ones(size(X,1)) v[2]*ones(size(X,1))]
bc = (:neu,:dir,:neu,:dir)

u  = X -> cos.(pi*X[:,1]).*cos.(pi*X[:,2])

ux = X -> -pi*sin.(pi*X[:,1]).*cos.(pi*X[:,2])
uy = X -> -pi*cos.(pi*X[:,1]).*sin.(pi*X[:,2])
sig = .001;

f  = X->  +2*sig*pi^2*u(X)+v[1]*ux(X)+v[2]*uy(X)

err = zeros(length(N))

bicg = (A,b; M=identity,tol=1e-10,maxIter=500,out=1)-> bicgstb(A,b,M1=identity,tol=tol,maxIter=maxIter,out=out)
times = zeros(length(N),3)
for j=1:2
	for k=1:length(N)
		Ainvs = (getJuliaSolver(), getMUMPSsolver(), ConvDiff.getBICGSTB(out=1))

	n      = [N[k];N[k]]
	M  = getRegularMesh(domain,n)
	xc = getCellCenteredGrid(M)

	pFor   = getConvDiffParam(M,V,sig=sig,gd=u,gn=ux,bc=bc,Ainv=Ainvs[j])
	
	dobs,pFor = getData(f(xc), pFor)
	
	utrue = u(xc)
	tic()
	err[k] = norm(pFor.Fields - utrue,Inf)/norm(utrue,Inf)
	times[k,j]= toq()
	println("err[$k] = $(err[k])")
end
@test all(abs.(diff(log2.(err))).>0.8)
end
