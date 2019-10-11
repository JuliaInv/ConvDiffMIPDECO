using jInv.LinearSolvers
using Test
using ConvDiffMIPDECO
using jInv.ForwardShare
using jInv.Mesh
using MUMPSjInv

domain = [0 1. 0 1]

N = [16;32;64;128]
v = randn(2)
v = v/norm(v)

u  = X -> cos.(pi*X[:,1]).*cos.(pi*X[:,2])
ux = X -> -pi*sin.(pi*X[:,1]).*cos.(pi*X[:,2])
uy = X -> -pi*cos.(pi*X[:,1]).*sin.(pi*X[:,2])
sig = .01;

f  = X ->  2*sig*pi^2*u(X)+v[1]*ux(X)+v[2]*uy(X)

err = zeros(length(N))

times = zeros(length(N),3)

for j=1:2
	for k=1:length(N)
		Ainvs = (getJuliaSolver(), getMUMPSsolver())


		nk      = [N[k];N[k]]
		Mk  = getRegularMesh(domain,nk)
		xck = getCellCenteredGrid(Mk)
		xnk = getNodalGrid(Mk)
		pFork = getConvDiffFEMParam(Mk,sig=sig,v=v,bc=(:dir,:neu,:neu,:neu),gd=u,Ainv=Ainvs[j])

		dobsk,pFork = getData(f(xck), pFork)
		# pFor.Fields -= mean(pFor.Fields)
		utruek = u(xnk)
		# utrue -= mean(utrue)
		times[k,j] = @elapsed begin
			err[k] = norm(pFork.Fields - utruek,Inf)/norm(utruek,Inf)
		end
		println(err[k])
end
@test all(abs.(diff(log2.(err))).>1.2)
end
