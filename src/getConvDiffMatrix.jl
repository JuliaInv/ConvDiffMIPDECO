export getConvDiffMatrix

"""
function getConvDiffMatrix(M::RegularMesh,sig::Float64,v,iddir::Array{Int,1},idneu::Array{Int,1}, 		  
	                      iddirn::Array{Int,1},idneun::Array{Int,1},idint::Array{Int,1})

builds convection diffusion matrix on interior nodes and returns contribution of dirichlet and neuman boundary.
"""
function getConvDiffMatrix(M::RegularMesh,sig::Float64,v,iddir::Array{Int,1},idneu::Array{Int,1}, iddirn::Array{Int,1},idneun::Array{Int,1},idint::Array{Int,1})
	
	# biuld matrix on padded mesh
	Lap  = getLaplacianMatrix(M)
	Conv = getConvectionMatrix(M,v)
	Apad = sig*Lap + Conv
	
	# projection for gost point neighbors
	I          = speye(prod(M.n+2))[:,idint]
	I[iddir,:] = -I[iddirn,:]  # dirichlect boundary points
	I[idneu,:] =  I[idneun,:]  # neuman boundary points
	
	# matrix for interior nodes
	A    = Apad[idint,:]*I
	
	# boundary nodes
	Adir = Apad[idint,iddir] 
	Aneu = Apad[idint,idneu]
	
	return A, Adir, Aneu
end
