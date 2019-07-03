using jInv.LinearSolvers
using Base.Test
using ConvDiff
using jInv.ForwardShare
using jInv.Mesh
using PyPlot
using jInvVis

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
    	
    	
		n      = [N[k];N[k]]
		M  = getRegularMesh(domain,n)
		xc = getCellCenteredGrid(M)
		xn = getNodalGrid(M)
		pFor = getConvDiffFEMParam(M,sig=sig,v=v,bc=(:dir,:neu,:neu,:neu),gd=u,Ainv=Ainvs[j])
		
		dobs,pFor = getData(f(xc), pFor)
		# pFor.Fields -= mean(pFor.Fields)
		utrue = u(xn)
		# utrue -= mean(utrue)
		tic()
		err[k] = norm(pFor.Fields - utrue,Inf)/norm(utrue,Inf)
		times[k,j]= toq()
		println(err[k])
end
@test all(abs.(diff(log2.(err))).>1.2)
end
