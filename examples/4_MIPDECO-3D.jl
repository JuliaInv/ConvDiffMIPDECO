using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MUMPSjInv
using MAT
using LinearAlgebra


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
mref       = zeros(M.nc)
reg        = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)
# reg        = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-8)

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

roundings = (s->naiveRounding(s),s->massPreservingRounding(s), s->objGapRedRounding(s,pMis,pInv))
roundingsNames = ["naive","massPreserving","gapReduction"];
neighborhoods = (s->trues(size(s)), s->dilation(s,pInv.MInv,1), s->dilation(s,pInv.MInv,2))
neighborhoodsNames = ["all","dilation1", "dilation2"];

nrow = length(roundings)
ncol = length(neighborhoods)

results = zeros(m[1],m[2],m[3],length(roundings),length(neighborhoods))
init = zeros(m[1],m[2],m[3],length(roundings),length(neighborhoods))
His = zeros(pInv.maxIter,3,length(roundings),length(neighborhoods))
times = zeros(length(roundings),length(neighborhoods))

for k1=3:length(roundings)
    for k2=1:length(neighborhoods)
        println("\n --- mipdecoHeuristic starting with $(roundingsNames[k1]) using $(neighborhoodsNames[k2])--")
        (t0,nsolve0) = (pFor.Ainv.solveTime, pFor.Ainv.nSolve)
        roundTime = @elapsed begin
            src0 = roundings[k1](srcRelaxed)
        end
        println("\t\ttime for rounding: $roundTime")
        dt  = pFor.Ainv.solveTime - t0
        println("\t\ttime for PDEsolves: $dt")
        dn = pFor.Ainv.nSolve - nsolve0
        println("\t\tnumber of PDEsolves: $dn")
        init[:,:,:,k1,k2] = src0

        times[k1,k2] = @elapsed begin
            mcTR,DcTR,flagTR,his = mipdecoHeuristic(src0,pInv,pMis,getNeighborhood=neighborhoods[k2],out=0)
            His[:,:,k1,k2] = his;
        end
        nIter = findlast(His[:,1,k1,k2].>0)


        println("\t\t initial obj. :\t\t$(His[1,1,k1,k2]) ")
        println("\t\t final obj. :\t\t$(His[nIter,1,k1,k2]) ")
        println("\t\t relative :\t\t\t$(His[nIter,1,k1,k2]/His[1,1,k1,k2]) ")
        println("\t\t no. PDE solves:\t$(2*count(His[:,1,k1,k2].>0)) ")
        println("\t\t runtime: \t\t\t$(times[k1,k2])")
        results[:,:,:,k1,k2] = mcTR
    end
end
res["SourcesMIPDECO"]=results
res["InitMIPDECO"]=init
res["HisMIPDECO"] = His
res["timesMIPDECO"] = times

matwrite(resFile,res)
