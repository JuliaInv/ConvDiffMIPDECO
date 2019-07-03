using Test
using ConvDiffMIPDECO
using jInv.Mesh

domain = [0 1. 0 1 0 1]
n  = [7 8 9]
M2D  = getRegularMesh(domain[1:4],n[1:2])

gx = X -> X[:,1]
gy = X -> X[:,2]
gz = X -> X[:,3]

xa,ya = getCellCenteredAxes(M2D)
(gl,gb,gr,gt) = getBoundaryCondition(gx, M2D)
@test all(gl.==0.0)
@test all(gb.==xa)
@test all(gr.==1.0)
@test all(gt.==xa)
@test length(gl)==n[2]
@test length(gb)==n[1]
@test length(gr)==n[2]
@test length(gt)==n[1]


(gl,gb,gr,gt) = getBoundaryCondition(gy, M2D)
@test all(gl.==ya)
@test all(gb.==0.0)
@test all(gr.==ya)
@test all(gt.==1.0)


M3D   = getRegularMesh(domain,n)
xa,ya,za = getCellCenteredAxes(M3D)
(gx1,gy1,gz1,gx2,gy2,gz2) = getBoundaryCondition(gx, M3D)
@test all(gx1.==0.0)
@test all(gx2.==1.0)
@test all(gy1.==repmat(xa,1,n[3]))
@test all(gy2.==repmat(xa,1,n[3]))
@test all(gz1.==repmat(xa,1,n[2]))
@test all(gz2.==repmat(xa,1,n[2]))


