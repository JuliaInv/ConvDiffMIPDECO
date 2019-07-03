using ConvDiff
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers


function runInv()
domain = [0 3. -0.5 1.5]
n      = [60 40]
M      = getRegularMesh(domain,n)
display(M)

ftrue = getFEMsource2D(M)
#pFor  = getConvDiffFEMParam(M,sig=.01,Ainv = getMUMPSsolver())
 pFor  = getConvDiffFEMParam(M,sig=.01,Ainv = getJuliaSolver())
dobs, = getData(vec(ftrue),pFor)
dobs += 0.4*randn(size(dobs))*mean(abs.(vec(dobs))) # add noise

sigback     = 0.0
Wt          = ones(prod(n+1))
pMis        = getMisfitParam(pFor,Wt,dobs,SSDFun);


# configure regularization
alpha	   	= 1e-2;
reg         = (m,mr,M,I=1.0) -> wTVReg(m,mr,M,eps=1e-6)
mref        = zeros(M.nc);
# configure optimization
HesPrec     = getSSORRegularizationPreconditioner(1.0,1e-15,50);
cgit       	= 5; 
pcgTol     	= 1e-1;
maxIter    	= 5;
minUpdate 	= 1e-3;

# use identity model
modFun      = identityMod
boundsLow  	= 0*ones(prod(n));
boundsHigh 	= 1*ones(prod(n));

# # use level set approach
modFun      = boundMod
boundsLow  	= -Inf*ones(prod(n));
boundsHigh 	= Inf *ones(prod(n));

maxStep		= 0.1*maximum(boundsHigh);

pInv = getInverseParam(M,modFun,reg,alpha,mref,
         boundsLow,boundsHigh,maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
         minUpdate=minUpdate, maxIter = maxIter,HesPrec=HesPrec);

mc,Dc,flag,His = projGNCG(0.0+0*ftrue[:],pInv,pMis);
fc = pInv.modelfun(mc)[1];

println("time grad misfit: $(sum(His.timeGradMisfit))")
println("time GN solve:    $(sum(His.timeLinSol))")
println("time misfit       $(sum(His.timeMisfit))")
println("time reg          $(sum(His.timeReg))")
end

runInv()

# using ProfileView
# Profile.clear()
# 
# @profile runInv()
# 
# ProfileView.view()
