module ConvDiffMIPDECO

using jInv.Mesh
using jInv.ForwardShare
using jInv.Utils
using jInv.LinearSolvers
using jInv.InverseSolve
using KrylovMethods
using LinearAlgebra
using SparseArrays
using Printf
using DSP

function getBICGSTB(;PC=:jac,maxIter=1000,out=0,tol=1e-10)
	bicg = (A,b; M=identity,tol=1e-10,maxIter=500,out=1)->
	bicgstb(A,b,M1=identity,tol=tol,maxIter=maxIter,out=out,tolRho=1e-60)
	return getIterativeSolver(bicg,PC=PC,maxIter=maxIter,out=out,tol=tol)
end

# files are organized in the src (containing functions) and test (containing some unit tests)
import jInv.ForwardShare.ForwardProbType

export ConvDiffParam, getConvDiffParam

"""
type ConvDiffParam <: ForwardProbType

description of stationary convection diffusion forward problem

(- sig*Laplacian + v dot GRAD) u = f
with boundary conditions du/dn = gn on Omega1, u=gd on Omega2

Construct an instance using getConvDiffParam(M,v,kwargs...)

Fields:
	M      - describes computational mesh
	v      - velocities at cell-centers
	sig    - scalar, >0
	bc     - vector describing boundary conditions
	Fields - stores fields, i.e., u
	Ainv   - factorization of PDE

"""
mutable struct ConvDiffParam <: ForwardProbType
	M :: RegularMesh
	A :: SparseMatrixCSC{Float64}
	P :: AbstractArray{Float64}
	sig::Float64
	bc::Array{Float64}
	Fields::Array{Float64,1}
	Ainv::AbstractSolver
end

"""
function getConvDiffParam(M,v)

constructs and returns ConvDiffParam

The PDE is discretized and factorized and stored in field Ainv.

Inputs:
	M            - mesh

Keyword arguments:
	v            - velocity, (vector or function)
	gd::Function - Dirichlet boundary conditions
	gn::Function - Neuman boundary condition
	P            - measurement operator
	sig          - viscosity
	bc           - description of boundary conditions
	Fields       - storing PDE solutions
	Ainv         - description of linear solver

"""
function getConvDiffParam(M::RegularMesh,v;
	gd::Function=X->zeros(size(X,1)), gn::Function=X->zeros(size(X,1)), P=Diagonal(ones(M.nc)),sig::Number=0.01,bc=(:dir,:neu,:neu,:neu),Fields=zeros(0),Ainv=getJuliaSolver())

	# get boundary conditions
	iddir, idneu, iddirn, idneun, idint = getBoundaryIndices(M,bc)

	A, Adir, Aneu = getConvDiffMatrix(M,sig,v,iddir,idneu,iddirn,idneun,idint)

	gdir, gneu = getBoundaryFuncVal(gd,gn,M,bc)

	bc = -2*Adir*gdir + Aneu*gneu

	return ConvDiffParam(M,A,P,sig,bc,Fields,Ainv)
end


include("getBoundaryCondition.jl")
include("getGhostIndices.jl")
include("getConvDiffMatrix.jl")
include("getBoundaryIndices.jl")
include("getBoundaryFuncVal.jl")
include("getDiffOps.jl")
include("getData.jl")
include("utils.jl")
include("ConvDiffFEMParam.jl")
include("FEM.jl")
include("getDataFEM.jl")
include("getConvDiffFEMConstraintsAMPL.jl")
include("mipdecoHeuristic.jl")
include("rounding.jl")
include("dilation.jl")
include("regularizers.jl")
end
