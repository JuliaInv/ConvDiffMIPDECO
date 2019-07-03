using ConvDiff
using Test
using jInv.Mesh
using MAT
using jInvVis
using PyPlot


domain = [0. 3. 0 1.]
n      = [60 20]-1
M      = getRegularMesh(domain,n)
@time Mass, Mass_const, SM = getFEMMatrices2D(M)

# matfile = matread("testFEMmatlab.mat")
# SMm         = matfile["SM"]
# Massm       = matfile["Mass"]
# Mass_constm = matfile["Mass_const"]

# @test all(SMm .== SM)
# @test all(Massm .== Mass)
# @test all(Mass_constm .== Mass_constm)

f = getFEMsource2D(M)
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
