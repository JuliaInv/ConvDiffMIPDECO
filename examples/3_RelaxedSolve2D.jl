using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MUMPSjInv
using LinearAlgebra
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

resFile = joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples",
                    "$(dataset)-$(reg)-noise-$(noiseLevel)-$(m[1])x$(m[2]).mat")
res = matread(resFile)
dobs = res["dobs"]
idx = 11;
alph = res["alphas"][idx]
println("solve using alpha = $(alph)")

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
# reg        = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)

## Configure optimization
maxIter    = 100
minUpdate  = 1e-3
HesPrec    = getSSORRegularizationPreconditioner(1.0,1e-15,50)
cgit       = 5
pcgTol     = 1e-1
modFun     = identityMod
boundsLow  = 0*ones(M.nc)
boundsHigh = 1*ones(M.nc)
maxStep	   = 0.1*maximum(boundsHigh)

## store the configuration
pInv       = getInverseParam(M,modFun,reg,alph,mref,
                             boundsLow,boundsHigh,maxStep=maxStep,
                            pcgMaxIter=cgit,pcgTol=pcgTol,minUpdate=minUpdate,maxIter=maxIter,
                            HesPrec=HesPrec);

mc = mref.+0.1

println("alpha:\t\t\t\t$(alph)")
println("no. cells:\t\t\t$(length(mc))")
println("no. nodes:\t\t\t$(size(P,2))")
println("no. measurements:\t$(size(P,1))")
runtimeRelaxed = @elapsed begin
    mc,Dc,flag,HisRelaxed = projGNCG(mc,pInv,pMis)
    uc = pFor.Fields;          # u variables
    wr = pInv.modelfun(mc)[1]; # w variables
end
println("runtime:\t$runtimeRelaxed")
# get number of PDE solves
nCGIter = cgit * length(HisRelaxed.hisLinSol)
nGNIter = count(HisRelaxed.Jc.>0)
println("no. PDE solves:\t$(pFor.Ainv.nSolve)")
println("time PDE factor:\t$(pFor.Ainv.facTime)")
println("time PDE solves:\t$(pFor.Ainv.solveTime)")


FieldsRelaxed        = uc
SourcesRelaxed       = wr
PredictedDataRelaxed = Dc

println("err(Fields) = $((norm(FieldsRelaxed - res["Fields"][:,idx]))/norm(res["Fields"][:,idx]))")
println("err(Sources) = $((norm(SourcesRelaxed - res["Sources"][:,idx]))/norm(res["Sources"][:,idx]))")
println("err(PredictedData) = $((norm(PredictedDataRelaxed - res["PredictedData"][:,idx]))/norm(res["PredictedData"][:,idx]))")

res["FieldsRelaxed"] = FieldsRelaxed
res["SourcesRelaxed"] = SourcesRelaxed
res["PredictedDataRelaxed"] = PredictedDataRelaxed
res["runtimeRelaxed"] = runtimeRelaxed
res["alphaRelaxed"] = alph
res["HisRelaxed"] = HisRelaxed

matwrite(resFile,res)
