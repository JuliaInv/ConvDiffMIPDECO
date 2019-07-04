using Test
using ConvDiffMIPDECO
using jInv.ForwardShare
using jInv.LinearSolvers
using KrylovMethods
using jInv.Mesh
using SparseArrays
using LinearAlgebra

domain = [0 1. 0 3. 0 2.]

N = [8;16;32]
v = randn(3,1)
v /= norm(v)
V = X -> ones(size(X,1))*v'
bc = (:neu,:dir,:dir,:neu,:dir,:dir)

u  = X -> cos.(pi*X[:,1]).*cos.(pi*X[:,2]).*cos.(pi*X[:,3])
ux = X -> -pi*sin.(pi*X[:,1]).*cos.(pi*X[:,2]).*cos.(pi*X[:,3])
uy = X -> -pi*cos.(pi*X[:,1]).*sin.(pi*X[:,2]).*cos.(pi*X[:,3])
uz = X -> -pi*cos.(pi*X[:,1]).*cos.(pi*X[:,2]).*sin.(pi*X[:,3])
sig = .05;

f  = X->  3*sig*pi^2*u(X)+v[1]*ux(X)+v[2]*uy(X)+v[3]*uz(X)

err = zeros(length(N))

times = zeros(length(N),3)
for j=1:3
	for k=1:length(N)
		Ainvs = (getJuliaSolver(), getMUMPSsolver(),ConvDiffMIPDECO.getBICGSTB())
		
		nk      = [N[k];N[k];N[k]]
		Mk  = getRegularMesh(domain,nk)
		xck = getCellCenteredGrid(Mk)
	
		pFork   = getConvDiffParam(Mk,V,sig=sig,gd=u,bc=bc,Ainv=Ainvs[j])
		
		dobsk,pFork = getData(vec(f(xck)), pFork)
		
		utruek = u(xck)
		times[k,j] = @elapsed begin
			err[k] = norm(pFork.Fields - utruek,Inf)/norm(utruek,Inf)
		end
		println("err[$k] = $(err[k])")
	end
	@test all(abs.(diff(log2.(err))).>0.7)
end