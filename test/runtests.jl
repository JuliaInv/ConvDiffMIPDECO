using LinearAlgebra
using SparseArrays


println("\t=== test MIPDECO julia package ===")
include("testGetBoundaryCondition.jl")
include("testGetBoundaryFuncVal.jl")
include("testGetGhostIndices.jl")
include("testConvDiffParam.jl")
include("testFEM.jl")
include("testFEM3D.jl")
include("testRun2D.jl")
include("testRegularizers.jl")
include("testConvDiff2D.jl")
include("testConvDiff3D.jl")
include("testFEM2Dmms.jl")
include("testFEM3Dmms.jl")
include("testGetConvDiffFEMConstraintsAMPL.jl")
println("\t=== all tests passed ===")
