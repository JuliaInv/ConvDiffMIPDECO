using Test
using LinearAlgebra
using SparseArrays
using Statistics


@testset "ConvDiffMIPDECO" begin
@testset "getBoundaryCondition" begin
	include("testGetBoundaryCondition.jl")
end
@testset "getBoundaryFuncVal" begin
	include("testGetBoundaryFuncVal.jl")
end
@testset "getGhostIndices" begin
	include("testGetGhostIndices.jl")
end
@testset "ConvDiffParam" begin
	include("testConvDiffParam.jl")
end
@testset "FEM" begin
	include("testFEM.jl")
end
@testset "FEM3D" begin
	include("testFEM3D.jl")
end
@testset "run2D" begin
	include("testRun2D.jl")
end
@testset "convDiff2D" begin
	include("testConvDiff2D.jl")
end
@testset "convDiff3D" begin
	include("testConvDiff3D.jl")
end
@testset "FEM2Dmms" begin
	include("testFEM2Dmms.jl")
end
@testset "FEM3Dmms" begin
	include("testFEM3Dmms.jl")
end
@testset "getConvDiffFEMConstraintsAMPL" begin
	include("testGetConvDiffFEMConstraintsAMPL.jl")
end
end
