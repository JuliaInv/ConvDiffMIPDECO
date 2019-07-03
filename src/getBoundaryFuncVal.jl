export getBoundaryFuncVal

"""
gdir, gneu = function getBoundaryFuncVal(gd::Function, gn::Function, M::RegularMesh,bc=(:dir,:neu,:neu,:neu))

evaluates functions gd and gn on boundary faces. Function input arguments must be arrays that store points rowwise, e.g.,

gd = X -> X[:,1].*X[:,2].*exp(X[:,3])
"""
function getBoundaryFuncVal(gd::Function, gn::Function, M::RegularMesh,bc=(:dir,:neu,:neu,:neu))
	# get boundary conditions
	gdx  = getBoundaryCondition(gd, M)
	gnx  = getBoundaryCondition(gn, M)
	gdir = zeros(0)
	gneu = zeros(0)
	for k=1:length(bc)
		if bc[k]==:dir
			gdir = [gdir; vec(gdx[k])]
		elseif bc[k]==:neu
			inout = (k<=M.dim) ? -1 : 1
			comp = mod.((1:2*M.dim) .- 1,M.dim) .+ 1
			gneu = [gneu; inout*M.h[comp[k]]*vec(gnx[k])]
		end
	end
	return gdir, gneu
end