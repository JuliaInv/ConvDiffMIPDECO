
using ConvDiffMIPDECO
using PyPlot
using jInv.Mesh
using jInv.ForwardShare
using jInv.InverseSolve
using jInvVis
using jInv.LinearSolvers
using MAT
using MUMPSjInv

ENV["DYLD_LIBRARY_PATH"]="/usr/lib/"

# filename= "2DmodelLShaped.mat"
filename= joinpath(dirname(pathof(ConvDiffMIPDECO)),"examples","Peaks2D.mat")
file = matread(filename)
wtrue = vec(file["W"]) ## Sources on finest mesh
domain = file["domain"]
m      = file["m"]

M      = getRegularMesh(domain,m)
viewImage2D(wtrue,M)

rec = rand(200,2) .*[domain[2] domain[4]]
plot(rec[:,1],rec[:,2],".r")


# build linear interpolation matrix from nodes to receiver locations
x1,x2 = getNodalAxes(M)
P = interpmat(x1,x2, rec);

# build param that holds forward problem
pFor = getConvDiffFEMParam(M,Ainv=getMUMPSsolver(),P=P);

dtrue, pFor = getData(vec(wtrue),pFor);
utrue = pFor.Fields

file["utrue"]  = reshape(utrue,M.n[1]+1,M.n[2]+1);
file["dtrue"]  = dtrue
file["rec"]    = rec
file["sig"]    = pFor.sig
file["v"]      = pFor.v

matwrite(filename,file)

viewImage2D(utrue,getPaddedMesh(M))
plot(rec[:,1],rec[:,2],".r")
title("true fields")
