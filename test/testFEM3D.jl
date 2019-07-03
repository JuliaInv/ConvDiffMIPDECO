using ConvDiff
using Base.Test
using jInv.Mesh
using MAT
using jInvVis
using PyPlot


domain = [0. 3. 0 1. 0 2.]
n      = 3*[7 9 12]-1
M      = getRegularMesh(domain,n)
Mass, Mass_const, SM = getFEMMatrices3D(M)

e = ones(prod(M.n+1))
@test abs(prod((domain[2:2:end]-domain[1:2:end])) - dot(e,Mass*e))/dot(e,Mass*e) < 1e-2

f = getFEMsource3D(M)
v = Mass_const*f
u = SM\v


# figure(1)
# subplot(2,3,1)
#  viewImage2D(f,M)
# title("source")
# colorbar()
# 
# subplot(2,3,2)
#  viewImage2D(matfile["rhs"],M)
# title("MATLAB source")
# colorbar()
# 
# subplot(2,3,3)
#  viewImage2D(f-matfile["rhs"],M)
# title("difference")
# colorbar()
# 
# subplot(2,3,4)
# Mn = getPaddedMesh(M)
# viewImage2D(u,Mn)
# colorbar()
# title("fields")
# 
# subplot(2,3,5)
# viewImage2D(matfile["u"],Mn)
# colorbar()
# title("MATLAB fields")
# 
# subplot(2,3,6)
# viewImage2D(matfile["u"]-u,Mn)
# colorbar()
# title("difference")
