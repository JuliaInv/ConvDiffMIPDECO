using Base.Test
using ConvDiff
using jInv.Mesh

domain = [0. 1. 0 1 0 1]
n      = [6 4 8]
M2D    = getRegularMesh(domain[1:4],n[1:2])

gx     = X -> X[:,1]
gy     = X -> ones(size(X,1))

xa,ya = getCellCenteredAxes(M2D)
gdir,gneu = getBoundaryFuncVal(gx, gy, M2D,(:dir,:dir,:dir,:dir))
@test isempty(gneu)
@test length(gdir) == sum(n[1:2])*2
@test all(gdir[1:n[2]].==0.0)
@test all(gdir[n[2]+(1:n[1])] .== xa)
@test all(gdir[sum(n[1:2])+(1:n[2])].==1.0)
@test all(gdir[n[1]+2*n[2]+(1:n[1])] .== xa)

gdir,gneu = getBoundaryFuncVal(gx, gy, M2D,(:neu,:neu,:neu,:neu))
@test isempty(gdir)
@test length(gneu) == sum(n[1:2])*2
@test all(gneu[1:n[2]].== -M2D.h[1])
@test all(gneu[n[2]+(1:n[1])] .== -M2D.h[2])
@test all(gneu[sum(n[1:2])+(1:n[2])].==M2D.h[1])
@test all(gneu[n[1]+2*n[2]+(1:n[1])] .== M2D.h[2])

M3D = getRegularMesh(domain,n)
xa,ya,za = getCellCenteredAxes(M3D)
gdir,gneu = getBoundaryFuncVal(gx, gy, M3D,(:dir,:dir,:dir,:dir,:dir,:dir))
@test isempty(gneu)
@test all(gdir[1:prod(n[2:3])].==0.0)
@test length(gdir) == 2*prod(n[2:3]) + 2*prod(n[[1,3]]) + 2*prod(n[1:2])
