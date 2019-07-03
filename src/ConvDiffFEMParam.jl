export ConvDiffFEMParam, getConvDiffFEMParam

"""
type ConvDiffFEMParam <: ForwardProbType

description of stationary convection diffusion forward problem

(- sig*Laplacian + v dot GRAD) u = f
with boundary conditions du/dn = 0 on Omega1, u=gd on Omega2

Here: Use a FEM discretization.

Construct an instance using getConvDiffFEMParam(M,kwargs...)

Fields:
	M         - describes computational mesh
	A         - stiffness matrix
	MassConst - mass matrix for piecewise constant rhs
	P         - receiver matrix
	sig       - scalar, >0 
	v         - velocity (assumed to be spatially homogeneous)
	ubc       - field from Dirichlet boundary conditions
	rbc       - right hand side contribution from Dirichlet bc
	iddir     - indices of nodes at which we have Dirichlet BCs
	idint     - unconstrained nodes
	Fields    - stores fields, i.e., u
	Ainv      - factorization of PDE
"""
type ConvDiffFEMParam <: ForwardProbType
	M :: RegularMesh
	A :: SparseMatrixCSC{Float64} # stiffness matrix
	MassConst::SparseMatrixCSC{Float64} # mass matrix for piecewise constant source
	P :: AbstractArray{Float64} # receiver matrix
	sig ::Float64        # viscosity
	v::Array{Float64}    # velocity
	ubc::Array{Float64}  # field from Dirichlet boundary condition
	rbc::Array{Float64}  # right hand side from Dirichlet boundary condition
	iddir::Array         # indices of nodes on Dirichlet boundary
	idint::Array         # indices of nodes not on Dirichlet boudnary
	Fields::Array{Float64,1}
	Ainv::AbstractSolver
end

"""
function getConvDiffFEMParam

constructor for ConvDiffFEMParam

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
"""
function getConvDiffFEMParam(M::RegularMesh;   v=eye(M.dim,1), gd::Function=X->zeros(size(X,1)), 
	                                           sig::Number=0.01,bc=(:dir,:neu,:neu,:neu,:neu,:neu), P=Diagonal(ones(prod(M.n+1))),
											   Fields=zeros(0),Ainv=getJuliaSolver())

	if M.dim==2
		Mass, Mass_const, A = getFEMMatrices2D(M,sig=sig,v=v)
	else
		Mass, Mass_const, A = getFEMMatrices3D(M,sig=sig,v=v)
	end
	
	iddir = getBoundaryIndicesFEM(M, bc)  # Dirichlet nodes
	idint = setdiff(1:size(A,1),iddir) # interior nodes
	Abc = A[idint,iddir]
	A = A[idint,idint]
	Mass_const = Mass_const[idint,:]
	
	xd         = getNodalGrid(M)[iddir,:]	
	ud         = gd(xd)
	ubc        = zeros(prod(M.n+1));
	ubc[iddir] = ud
	rbc        = Abc*ud;
	
	return ConvDiffFEMParam(M,A,Mass_const,P,sig,v,ubc,rbc,iddir,idint,Fields,Ainv)
end
