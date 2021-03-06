using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MUMPSjInv
using MAT
using LinearAlgebra
using SparseArrays

# filename= "2DmodelLShaped.mat"
filename= joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples","Sources3D")
file = matread(filename * ".mat")

# problem description
domain = file["domain"]
mf     = round.(Int64,file["m"])
v      = file["v"]
sig    = file["sig"]

# data and receiver location
dtrue  = file["dtrue"]
rec    = file["rec"]


Mfine = getRegularMesh(domain,mf);
x1c,x2c,x3c = getCellCenteredAxes(Mfine)
rec3D = [kron(ones(length(x3c)),rec) kron(x3c,ones(size(rec,1)))]



# set up experiments
Meshes       = [ getRegularMesh(domain,[96 48 48]);];
alphas       = exp10.(range(0,stop=-6,length=30))
regularizersName = ["wTVReg"];
regularizers = [wTVReg];
noiseLevel = 0.1;

for mk = 1:length(Meshes)
	M      = Meshes[mk]
	x1,x2,x3 = getNodalAxes(M)
	P = interpmat(x1,x2,x3, rec3D);
	for rk = 1:length(regularizers)
		regFun= regularizers[rk]
		regName= regularizersName[rk]

		# build linear interpolation matrix from nodes to receiver locations

		dnoise = randn(size(dtrue));
		dnoise /= norm(dnoise)
		dnoise *= noiseLevel*norm(dtrue)
		dobs = dtrue + dnoise

		pFor = getConvDiffFEMParam(M,v=v,sig=sig,P=P,Ainv=getMUMPSsolver());

		## configure misfit
		Wt         = ones(size(dobs))/sqrt(mf[3])
		sigback    = 0.0
		pMis       = getMisfitParam(pFor,Wt,dobs,SSDFun)

		## Configure regularization
		mref       = zeros(M.nc)
		alpha      = alphas[1]
		reg        = (m,mr,M,I=1.0) -> regFun(m,mr,M,eps=1e-8)
		# reg        = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)

		## Configure optimization
		maxIter    = 20
		minUpdate  = 1e-3
		HesPrec    = getSSORRegularizationPreconditioner(1.0,1e-15,50)
		cgit       = 5
		pcgTol     = 1e-1
		modFun     = identityMod
		boundsLow  = 0*ones(M.nc)
		boundsHigh = 1*ones(M.nc)
		maxStep	   = 0.1*maximum(boundsHigh)

		## store the configuration
		pInv       = getInverseParam(M,modFun,reg,alpha,mref,
		                             boundsLow,boundsHigh,maxStep=maxStep,
		                            pcgMaxIter=cgit,pcgTol=pcgTol,minUpdate=minUpdate,maxIter=maxIter,
		                            HesPrec=HesPrec);

		mc = mref.+0.1

		Mis    = zeros(length(alphas))
		Reg    = zeros(length(alphas))
		Fields = zeros(size(P,2),length(alphas))
		Sources = zeros(length(mc),length(alphas))
	    PredictedData   = zeros(size(P,1),length(alphas))

		nPDEsolves = 0;

		runtimeLcurve = @elapsed begin
			for ak = 1:length(alphas)
				pInv.alpha = alphas[ak]
				println("\n\n--- Inversion on $(M.n) grid with $(regName) and alpha=$(pInv.alpha) --- ")
				# mc = max.(.01,min.(.99,mc))
				mc,Dc,flag,His = projGNCG(mc,pInv,pMis)
				uc = pFor.Fields;          # u variables
				wr = pInv.modelfun(mc)[1]; # w variables

				Fields[:,ak] = uc;
				Sources[:,ak] = wr
				PredictedData[:,ak] = Dc

				iter = findlast(His.F .> 0)
				Mis[ak] = His.F[iter]
				Reg[ak] = His.Rc[iter]/alphas[ak]
		    end
		end
		println("no. PDE solves:\t$(pFor.Ainv.nSolve)")
		println("time PDE factor:\t$(pFor.Ainv.facTime)")
		println("time PDE solves:\t$(pFor.Ainv.solveTime)")
		println("\n overall runtime: $runtimeLcurve")

		results = Dict("n"=>M.n, "Mis"=>Mis, "Reg"=>Reg, "alphas"=> alphas,
		"noiseLevel"=>noiseLevel,"Fields"=> Fields, "Sources" => Sources,
		"PredictedData" => PredictedData, "dobs" => dobs,"runtime"=>runtimeLcurve)
		resFile = "$(filename)-$regName-noise-$(noiseLevel)-$(M.n[1])x$(M.n[2])x$(M.n[3]).mat"
		matwrite(resFile,results)
	end
end
