using Test
using ConvDiffMIPDECO
using jInv.Mesh
using LinearAlgebra

domain = [0 1. 0 1 0 1]
n  = [7 8 9]
M2D  = getRegularMesh(domain[1:4],n[1:2])

idg,idgn,idint = getGhostIndices(n[1:2])

xn = getNodalGrid(getPaddedMesh(M2D))

@test norm(xn[idg[1],1].- (-M2D.h[1]/2))<1e-15
@test norm(xn[idg[2],2].- (-M2D.h[2]/2))<1e-15
@test norm(xn[idg[3],1].- (1+M2D.h[1]/2))<1e-15
@test norm(xn[idg[4],2].- (1+M2D.h[2]/2))<1e-15
@test norm(xn[idgn[1],1].- (M2D.h[1]/2))<1e-15
@test norm(xn[idgn[2],2].- (M2D.h[2]/2))<1e-15
@test norm(xn[idgn[3],1].- (1-M2D.h[1]/2))<1e-15
@test norm(xn[idgn[4],2].- (1-M2D.h[2]/2))<1e-15
@test length(idgn[1])==n[2]
@test length(idgn[2])==n[1]
@test length(idgn[3])==n[2]
@test length(idgn[4])==n[1]
@test length(idg[1])==n[2]
@test length(idg[2])==n[1]
@test length(idg[3])==n[2]
@test length(idg[4])==n[1]
@test length(idint)==prod(n[1:2])


M3D  = getRegularMesh(domain,n)

idg,idgn,idint = getGhostIndices(n)

xn = getNodalGrid(getPaddedMesh(M3D))

@test norm(xn[idg[1],1].- (-M3D.h[1]/2))  < 1e-15
@test norm(xn[idg[2],2].- (-M3D.h[2]/2))  < 1e-15
@test norm(xn[idg[3],3].- (-M3D.h[3]/2))  < 1e-15
@test norm(xn[idg[4],1].- (1+M3D.h[1]/2)) < 1e-15
@test norm(xn[idg[5],2].- (1+M3D.h[2]/2)) < 1e-15
@test norm(xn[idg[6],3].- (1+M3D.h[3]/2)) < 1e-15
@test norm(xn[idgn[1],1].- (M3D.h[1]/2))  < 1e-15
@test norm(xn[idgn[2],2].- (M3D.h[2]/2))  < 1e-15
@test norm(xn[idgn[3],3].- (M3D.h[3]/2))  < 1e-15
@test norm(xn[idgn[4],1].- (1-M3D.h[1]/2))  < 1e-15
@test norm(xn[idgn[5],2].- (1-M3D.h[2]/2))  < 1e-15
@test norm(xn[idgn[6],3].- (1-M3D.h[3]/2))  < 1e-15

@test length(idgn[1])==prod(n[2:3])
@test length(idgn[2])==prod(n[[1,3]])
@test length(idgn[3])==prod(n[1:2])
@test length(idgn[4])==prod(n[2:3])
@test length(idgn[5])==prod(n[[1,3]])
@test length(idgn[6])==prod(n[1:2])

@test length(idg[1])==prod(n[2:3])
@test length(idg[2])==prod(n[[1,3]])
@test length(idg[3])==prod(n[1:2])
@test length(idg[4])==prod(n[2:3])
@test length(idg[5])==prod(n[[1,3]])
@test length(idg[6])==prod(n[1:2])

@test length(idint)==prod(n)



