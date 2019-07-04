
using ConvDiffMIPDECO
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using MAT


filename= joinpath(dirname(pathof(ConvDiffMIPDECO)),"..","examples","Sources3D.mat")
file = matread(filename)
wtrue = vec(file["W"]) ## Sources on finest mesh
domain = file["domain"]
m      = vec(file["m"])

M      = getRegularMesh(domain,m)

rec = rand(100,2) .*[domain[2] domain[4] ]

Ainv = getMUMPSsolver()

# build linear interpolation matrix from nodes to receiver locations
x1c,x2c,x3c = getCellCenteredAxes(M)
rec3D = [kron(ones(length(x3c)),rec) kron(x3c,ones(size(rec,1)))]

x1,x2,x3 = getNodalAxes(M)
P = interpmat(x1,x2,x3, rec3D);

# build param that holds forward problem
pFor = getConvDiffFEMParam(M,Ainv=Ainv,P=P);
dtrue, pFor = getData(vec(wtrue),pFor);
utrue = pFor.Fields


file["utrue"]  = reshape(utrue,tuple(M.n.+1...));
file["dtrue"]  = dtrue
file["rec"]    = rec
file["sig"]    = pFor.sig
file["v"]      = pFor.v

matwrite(filename,file)
