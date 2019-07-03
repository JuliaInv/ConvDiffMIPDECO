using ConvDiff
using Base.Test
using jInv.Mesh
using jInvVis

println("=== 2D test === ")

domain = [0. 3. 0 1. ]
n      = 2*[7 9 ]
M      = getRegularMesh(domain,n)
bc = (:dir,:neu,:neu,:neu)
gd = X->sin.(X[:,2]*pi)
pFor = getConvDiffFEMParam(M,bc=bc,gd=gd)

A,B,b = getConvDiffFEMConstraintsAMPL(M,gd=gd)

mc = rand(M.nc); 
dobs = getData(mc,pFor)
utrue = pFor.Fields;
utest = A\(B*mc+b); # -->    A *u = B*m + b
@test norm(utrue-utest)/norm(utrue) < 1e-10


println("=== 3D test === ")

domain = [0. 3. 0 1. 0 2.]
n      = 2*[7 9 12]
M      = getRegularMesh(domain,n)
bc = (:dir,:neu,:neu,:neu,:neu,:neu)
gd = X->sin.(X[:,2]*pi)
pFor = getConvDiffFEMParam(M,bc=bc,gd=gd)

A,B,b = getConvDiffFEMConstraintsAMPL(M,gd=gd)

mc = rand(M.nc); 
dobs = getData(mc,pFor)
utrue = pFor.Fields;
utest = A\(B*mc+b);
@test norm(utrue-utest)/norm(utrue) < 1e-10
