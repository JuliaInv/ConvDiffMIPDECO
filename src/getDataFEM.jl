export getData
import jInv.ForwardShare.getData
function getData(m,pFor::ConvDiffFEMParam)
	pFor.Ainv.doClear = 0 # PDE does not depend on m, so we shouldn't clear i
	
	# solve PDE for non-Dirichlet nodes
	Ut,pFor.Ainv      = solveLinearSystem(pFor.A, pFor.MassConst*m - pFor.rbc, pFor.Ainv)

	U                 = pFor.ubc
	U[pFor.idint]     = Ut
	dobs 		      = pFor.P*U
	pFor.Fields       = U
	
	return dobs,pFor
end

import jInv.ForwardShare.getSensMatVec
function getSensMatVec(v::Vector{Float64},m::Vector{Float64},pFor::ConvDiffFEMParam)
	dUt = solveLinearSystem(pFor.A, pFor.MassConst*v, pFor.Ainv)[1]
	dU             = zeros(prod(pFor.M.n.+1),1)
	dU[pFor.idint] = dUt
	return  pFor.P*dU
end

import jInv.ForwardShare.getSensTMatVec
function getSensTMatVec(v::Vector{Float64},m::Vector{Float64},pFor::ConvDiffFEMParam)
	return  pFor.MassConst'*solveLinearSystem(pFor.A,(pFor.P'*v)[pFor.idint], pFor.Ainv,1)[1]
end
