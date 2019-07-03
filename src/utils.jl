export getPaddedMesh, jInvObj

"""
function getPaddedMesh(M::RegularMesh)

return regular mesh padded by h/2 on all sided and with 
one cell added, so that the cell-centered points on M
are nodal points on the padded mesh.
"""
function getPaddedMesh(M::RegularMesh)
	domainPad = copy(M.domain); 
	domainPad[1:2:end] -= M.h/2
	domainPad[2:2:end] += M.h/2
	return getRegularMesh(domainPad,M.n+1)
end



"""
function jInvObj(mc,pMis,pInv)

evaluates the jInv objective function defined by pMis and pInv
"""
function jInvObj(mc::BitArray{1},pMis,pInv)
   return jInvObj(Array{Float64}(mc),pMis,pInv); 
end
function jInvObj(mc::Array{Float64},pMis,pInv)
    
    sig,dsig              = pInv.modelfun(mc)
	tic()
	Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis,true)
	to = toq()
	timeF = to
	dF                    = dsig'*dF
	
	# Compute regularizer
	tic()
	R,dR,d2R = computeRegularizer(pInv.regularizer,mc,pInv.mref,pInv.MInv,pInv.alpha)
	to = toq()
	timeR = to
	
	# Objective function
	Jc = F  + R
	gc = dF + dR
    
    return Jc,gc,Dc,F,R,timeF,timeR
end