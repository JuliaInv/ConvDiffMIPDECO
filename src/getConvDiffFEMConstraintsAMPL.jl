export getConvDiffFEMConstraintsAMPL
"""
function getConvDiffFEMConstraintsAMPL

Builds matrices to model constraints to be used in AMPL. 
To this end, we take the FEM discretization on all nodes of 
the mesh and incorporate the Dirichlet boundary condition. 
This yields the constraints

Stiffness * u = Mass*w     (on non-Dirichlet nodes)
            u = gd         (on Dirichled nodes)

We represent these constraints as

A * u = B*w + b, where

A = [Stiffness; I[dir,:]], B = [Mass; zeros], b = [zeros; gd(dir)]

!! Note: The linear system arising from this approach will not be symmetric. 
!!       It can be made so by eliminating the Dirichlet nodes. This is done in 
!!       getConvDiffFEMParam, which is preferable in terms of efficiency

Required Input:
	M         - describes computational mesh

Keyword Arguments:
	sig       - scalar, >0 
	v         - velocity (assumed to be spatially homogeneous)
	gd        - function for Dirichlet boundary conditions
	bc        - Tuple, defining the boundary conditions
	P         - receiver matrix
	Fields    - stores fields, i.e., u
	Ainv      - factorization of PDE

Outputs:
    A = [Stiffness; I[dir,:]] - FEM discretization of constraints on non-Dirichlet nodes
    B = [Mass; zeros]         - mass matrix to model source contribution to non-Dirichlet nodes
    b = [0; gd(dir)]          - vector to encode Dirichlet conditions.

"""
function getConvDiffFEMConstraintsAMPL(M::RegularMesh; v=eye(M.dim,1), gd::Function=X->zeros(size(X,1)), 
	                                           sig::Number=0.01,bc=(:dir,:neu,:neu,:neu,:neu,:neu))

# get mass and stiffness matrices for all nodes
if M.dim==2
	Mass, Massc, Afull = getFEMMatrices2D(M,sig=sig,v=v)
else
	Mass, Massc, Afull = getFEMMatrices3D(M,sig=sig,v=v)
end

# get indices of nodes on Dirichlet boundary and indices of interior nodes
iddir = ConvDiff.getBoundaryIndicesFEM(M, bc)  # Dirichlet nodes
idint = setdiff(1:size(Afull,1),iddir) # interior nodes

# take constraints from interior nodes and add rows for Dirichlet condition. 
# note that rows correspond to test functions, which are not needed on Dirichlet 
# boundary
Id = speye(size(Afull,2))
#A  = [Afull[idint,:];  Id[iddir,:]]
#B  = [Massc[idint,:]; spzeros(length(iddir),size(Massc,2))]
A  = [Id[iddir,:];Afull[idint,:]]
B  = [spzeros(length(iddir),size(Massc,2)); Massc[idint,:]]
xn = getNodalGrid(M)
#b  = [zeros(length(idint)); gd(xn[iddir,:])]
b  = [ gd(xn[iddir,:]);zeros(length(idint))]
return A,B,b
end
