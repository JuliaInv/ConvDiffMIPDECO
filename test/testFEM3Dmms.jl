using Base.Test
using ConvDiff
using jInv.ForwardShare
using jInv.LinearSolvers
using KrylovMethods
using jInv.Mesh

domain = [0 1. 0 3. 0 2.]

N = [8;16;32;]
v = X -> [ones(size(X,1)) zeros(size(X,1)) zeros(size(X,1))]
bc = (:dir,:neu,:neu,:neu,:neu,:neu)

u  = X -> cos.(pi*X[:,1]).*cos.(pi*X[:,2]).*cos.(pi*X[:,3])
ux = X -> -pi*sin.(pi*X[:,1]).*cos.(pi*X[:,2]).*cos.(pi*X[:,3])
uy = X -> -pi*cos.(pi*X[:,1]).*sin.(pi*X[:,2]).*cos.(pi*X[:,3])
uz = X -> -pi*cos.(pi*X[:,1]).*cos.(pi*X[:,2]).*sin.(pi*X[:,3])
sig = .01;

f  = X->  3*sig*pi^2*u(X)+ux(X)

err = zeros(length(N))

times = zeros(length(N),3)
pFor=1;
for j=1:2
	for k=1:length(N)
		Ainvs = (getJuliaSolver(), getMUMPSsolver())
	
		n      = [N[k];N[k];N[k]]
		M  = getRegularMesh(domain,n)
		xc = getCellCenteredGrid(M)
		xn = getNodalGrid(M)
		pFor = getConvDiffFEMParam(M,sig=sig,gd=u,bc=bc,Ainv=Ainvs[j])
		
		dobs,pFor = getData(f(xc), pFor)
		pFor.Fields -= mean(pFor.Fields)
		utrue = u(xn)
		utrue -= mean(utrue)
		tic()
		err[k] = norm(pFor.Fields - utrue,Inf)/norm(utrue,Inf)
		times[k,j]= toq()
		
		println(err[k])
	end
	@test all(abs.(diff(log2.(err))).>1.4)
end