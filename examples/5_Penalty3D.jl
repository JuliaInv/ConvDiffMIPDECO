using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MUMPSjInv
using MAT

# read results from l-curve
dataset    = "Sources3D"
reg        = "wTVReg"
m          = [96 48 48]
noiseLevel = 0.1

filename= joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples","$(dataset).mat")
data = matread(filename)
domain = data["domain"]
v      = data["v"]
sig    = data["sig"]
dtrue  = data["dtrue"]
rec    = data["rec"];
mf     = round.(Int64,data["m"])

M = getRegularMesh(domain,m)

resFile = joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples",
                 "$(dataset)-$(reg)-noise-$(noiseLevel)-$(m[1])x$(m[2])x$(m[3]).mat")
res = matread(resFile)
dobs = res["dobs"]
srcRelaxed = res["SourcesRelaxed"]
alphaRelaxed = res["alphaRelaxed"]

resFilePenalty = joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples",
                 "$(dataset)-$(reg)-noise-$(noiseLevel)-$(m[1])x$(m[2])x$(m[3])-penalty.mat")


println("solve using alpha = $(alphaRelaxed)")

noiseLevel = norm(dobs-dtrue)/norm(dtrue)
println("estimated noise level: $noiseLevel")
if noiseLevel > 0.5
    error("noise level too large, double check that data was loaded correctly. ")
end

# build inverse problem
Mfine = getRegularMesh(domain,vec(data["m"]));
x1c,x2c,x3c = getCellCenteredAxes(Mfine)
rec3D = [kron(ones(length(x3c)),rec) kron(x3c,ones(size(rec,1)))]
# build linear interpolation matrix from nodes to receiver locations
x1,x2,x3 = getNodalAxes(M)
P = interpmat(x1,x2,x3, rec3D);
pFor = getConvDiffFEMParam(M,v=v,sig=sig,P=P,Ainv=getMUMPSsolver());
## configure misfit
Wt         = ones(size(dobs))/sqrt(mf[3])
sigback    = 0.0
pMis       = getMisfitParam(pFor,Wt,dobs,SSDFun)

## Configure regularization
mref       = zeros(M.nc,2)
regTV      = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)
regP       = (m,mr,M,I=1.0) -> mipReg(m,mr,M)
reg        = [regTV;regP]

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
pInv       = getInverseParam(M,modFun,reg,[alphaRelaxed;0.0],mref,
                             boundsLow,boundsHigh,maxStep=maxStep,
                            pcgMaxIter=cgit,pcgTol=pcgTol,minUpdate=minUpdate,maxIter=maxIter,
                            HesPrec=HesPrec);

mc = mref[:,1].+0.1
# viewImage2D(mc,M)
println("alpha:\t\t\t\t$(pInv.alpha)")
println("no. cells:\t\t\t$(length(mc))")
println("no. nodes:\t\t\t$(size(P,2))")
println("no. measurements:\t$(size(P,1))")
runtimeRelaxed = @elapsed begin
    mc,Dc,flag,HisRelaxed = projGNCG(mc,pInv,pMis)
    uc = pFor.Fields;          # u variables
    wr = pInv.modelfun(mc)[1]; # w variables
end

println("\n\n--Penalty method--\n")
Imax = 30
epsilon = 1e-4
betaMin = 1e-6
pInv.alpha[2] = betaMin
pInv.maxIter = 10
Ms = zeros(length(mc),Imax+1)
iter=1;
Ms[:,1] = mc;
HIS = Array{Any}(undef,Imax)
runtimePenalty = @elapsed begin
for iter=1:Imax
    mround,Dc,flag,HIS[iter] = projGNCG(Ms[:,iter],pInv,pMis)
    uc = pFor.Fields;          # u variables
    wr = pInv.modelfun(mround)[1]; # w variables

    Ms[:,iter+1] = mround;
    pInv.alpha[2] *= 2
end
end
println("runtime:\t$runtimePenalty")
println("no. PDE solves:\t$(pFor.Ainv.nSolve)")
println("time PDE factor:\t$(pFor.Ainv.facTime)")
println("time PDE solves:\t$(pFor.Ainv.solveTime)")


res["HIS"] = HIS
res["Ms"] = Ms
res["runtimeRelaxed"] = runtimeRelaxed
res["runtimePenalty"] = runtimeRelaxed
res["alphaPenalty"] = pInv.alpha
res["HisRelaxed"] = HisRelaxed

matwrite(resFilePenalty,res)
