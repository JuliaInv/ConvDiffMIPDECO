
using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MAT

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

resFile = joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples","$(dataset)-$(reg)-noise-$(noiseLevel)-$(m[1])x$(m[2]).mat")
res = matread(resFile)
dobs = res["dobs"]

idx = 11;
alpha = res["alphas"][idx]
println("solve using alpha = $(alpha)")

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
# reg        = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)

## Configure optimization
maxIter    = 10
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


roundings = (s->naiveRounding(s),s->massPreservingRounding(s), s->objGapRedRounding(s,pMis,pInv))
neighborhoods = (s->trues(size(s)), s->dilation(s,pInv.MInv,1), s->dilation(s,pInv.MInv,2))

nrow = length(roundings)
ncol = length(neighborhoods)
pInv.maxIter=50

results = zeros(m[1],m[2],length(roundings),length(neighborhoods))
init = zeros(m[1],m[2],length(roundings),length(neighborhoods))
His = zeros(pInv.maxIter,3,length(roundings),length(neighborhoods))
times = zeros(length(roundings),length(neighborhoods))
for k1=1:length(roundings)
    for k2=1:length(neighborhoods)
        src0 = roundings[k1](res["Sources"][:,idx])
        init[:,:,k1,k2] = src0
        times[k1,k2] = @elapsed begin
        mcTR,DcTR,flagTR,his = mipdecoHeuristic(src0,pInv,pMis,getNeighborhood=neighborhoods[k2],out=0)
        results[:,:,k1,k2] = mcTR
        His[:,:,k1,k2] = his;
		end

    end
end



res["MIPDECO"]=results
res["MIPDECOinit"]=init
res["His"] = His
res["times"] = times
matwrite(resFile,res)
