using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MUMPSjInv
using MAT
# using jInvVis

# read results from l-curve
dataset    = "Peaks2D"
reg        = "wTVReg"
m          = [256 128]
noiseLevel = 0.1

filename= joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples","Peaks2D.mat")
data = matread(filename)
domain = data["domain"]
v      = data["v"]
sig    = data["sig"]
dtrue  = data["dtrue"]
rec    = data["rec"];

M = getRegularMesh(domain,m)

resFile = joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples",
                    "$(dataset)-$(reg)-noise-$(noiseLevel)-$(m[1])x$(m[2])-penalty.mat")
res = matread(resFile)
dobs = res["dobs"]



srcRelaxed   = res["SourcesRelaxed"]
alphaRelaxed = res["alphaRelaxed"]
Ms           = res["Ms"]

println("solve using alpha = $(alphaRelaxed)")
noiseLevel = norm(dobs-dtrue)/norm(dtrue)
println("estimated noise level: $noiseLevel")
if noiseLevel > 0.5
    error("noise level too large, double check that data was loaded correctly. ")
end


# build inverse problem
x1,x2 = getNodalAxes(M)
P = interpmat(x1,x2, rec);
pFor = getConvDiffFEMParam(M,v=v,sig=sig,P=P,Ainv=getMUMPSsolver());

## configure misfit
Wt         = ones(size(dobs))
sigback    = 0.0
pMis       = getMisfitParam(pFor,Wt,dobs,SSDFun)

## Configure regularization
mref       = zeros(M.nc)
reg        = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)

## Configure optimization
maxIter    = 50
minUpdate  = 1e-3
HesPrec    = getSSORRegularizationPreconditioner(1.0,1e-15,50)
cgit       = 5
pcgTol     = 1e-1
modFun     = identityMod
boundsLow  = 0*ones(M.nc)
boundsHigh = 1*ones(M.nc)
maxStep	   = 0.1*maximum(boundsHigh)

## store the configuration
pInv       = getInverseParam(M,modFun,reg,alphaRelaxed,mref,
                             boundsLow,boundsHigh,maxStep=maxStep,
                            pcgMaxIter=cgit,pcgTol=pcgTol,minUpdate=minUpdate,maxIter=maxIter,
                            HesPrec=HesPrec);


neighborhoods = (s->trues(size(s)), s->dilation(s,pInv.MInv,1), s->dilation(s,pInv.MInv,1,0.1,0.9))
neighborhoodsNames = ["all","dilation1", "dilation2"];

nrow = size(Ms,2)
ncol = length(neighborhoods)

results = zeros(m[1],m[2],nrow,length(neighborhoods))
init = zeros(m[1],m[2],nrow,length(neighborhoods))
His = zeros(pInv.maxIter,3,nrow,length(neighborhoods))
times = zeros(nrow,length(neighborhoods))


println("solve using alpha = $(pInv.alpha)")

for k1=1:nrow
    for k2=1:length(neighborhoods)
        println("\n --- mipdecoHeuristic starting with Ms[:,$k1] using $(neighborhoodsNames[k2])--")
        (t0,nsolve0) = (pFor.Ainv.solveTime, pFor.Ainv.nSolve)

        src0 = naiveRounding(Ms[:,k1])
        init[:,:,k1,k2] = src0


        times[k1,k2] = @elapsed begin
            mcTR,DcTR,flagTR,his = mipdecoHeuristic(src0,pInv,pMis,getNeighborhood=neighborhoods[k2],out=0)
            results[:,:,k1,k2] = mcTR
            His[:,:,k1,k2] = his;
		end
        nIter = findlast(His[:,1,k1,k2].>0)

        # p1 = viewImage2D(Ms[:,k1],M)
        # p2 = viewImage2D(src0,M)
        # p3 = viewImage2D(mcTR,M)
        # p4 = viewImage2D((mcTr-src0),M)
        # pt = plot(p1,p2,p3,p4,layout=(4,1))
        # display(pt)

        println("\t\t initial obj. :\t\t$(His[1,1,k1,k2]) ")
        println("\t\t final obj. :\t\t$(His[nIter,1,k1,k2]) ")
        println("\t\t relative :\t\t\t$(His[nIter,1,k1,k2]/His[1,1,k1,k2]) ")
        println("\t\t no. PDE solves:\t$(2*nIter) ")
        println("\t\t runtime: \t\t\t$(times[k1,k2])")
    end
end

res["SourcesMIPDECOPenalty"]=results
res["InitMIPDECOPenalty"]=init
res["HisMIPDECOPenalty"] = His
res["timesMIPDECOPenalty"] = times
matwrite(resFile,res)
