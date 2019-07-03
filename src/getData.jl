
import jInv.ForwardShare.getData
function getData(m,pFor::ConvDiffParam)
	pFor.Ainv.doClear = 0 # PDE does not depend on m, so we shouldn't clear it

	U,pFor.Ainv    = solveLinearSystem(pFor.A, m + pFor.bc, pFor.Ainv)
	dobs = pFor.P*U
	pFor.Fields = U
	
	return dobs,pFor
end

import jInv.ForwardShare.getSensMatVec
function getSensMatVec(v::Vector{Float64},m::Vector{Float64},pFor::ConvDiffParam)
	return pFor.P*solveLinearSystem(pFor.A, v, pFor.Ainv)[1]
end

import jInv.ForwardShare.getSensTMatVec
function getSensTMatVec(v::Vector{Float64},m::Vector{Float64},pFor::ConvDiffParam)
	return solveLinearSystem(pFor.A, pFor.P'*v, pFor.Ainv,1)[1]
end
