#!/bin/sh

julia 3_Lcurve2D.jl > 3_Lcurve2D.txt
julia 4_RelaxedSolve2D.jl > 4_RelaxedSolve2D.txt
julia 5_MIPDECO2D.jl > 5_MIPDECO2D.txt

julia 3_Lcurve3D.jl > 3_Lcurve3D.txt
julia 4_RelaxedSolve3D.jl > 4_RelaxedSolve3D.txt
julia 5_MIPDECO3D.jl > 5_MIPDECO3D.txt
