#!/bin/sh

# matlab getPeaksModel2D.m
# matlab getModel3D.m

# julia 1_GenerateData.jl > 1_GenerateData.txt
# julia 2_Lcurve2D.jl > 2_Lcurve2D.txt
# julia 3_RelaxedSolve2D.jl > 3_RelaxedSolve2D.txt
# julia 4_MIPDECO2D.jl > 4_MIPDECO2D.txt
# julia 5_Penalty2D.jl > 5_Penalty2D.txt
# julia 6_MIPDECO2DPenalty.jl > 6_MIPDECO2DPenalty.txt


# julia 1_GenerateData.jl > 1_GenerateData.txt
# julia 2_Lcurve3D.jl > 2_Lcurve3D.txt
#julia 3_RelaxedSolve3D.jl > 3_RelaxedSolve3D.txt
julia 4_MIPDECO-3D.jl > 4_MIPDECO-3D.txt
julia 5_Penalty3D.jl > 5_Penalty3D.txt
julia 6_MIPDECO3DPenalty.jl > 6_MIPDECO3DPenalty.txt
