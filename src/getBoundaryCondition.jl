export getBoundaryCondition

"""
function getBoundaryCondition(gd::Function, M::AbstractMesh)

evaluate function gd on boundary faces of mesh M

"""
function getBoundaryCondition(gd::Function, M::AbstractMesh)
	if M.dim==2
		xf1,xf2 = getFaceGrids(M)
	
		gf1  = reshape(gd(xf1),tuple(M.n+[1,0]...))
		gdl  = gf1[1,:]
		gdr  = gf1[end,:]

		gf2  = reshape(gd(xf2),tuple(M.n+[0,1]...))
		gdb  = gf2[:,1]
		gdt  = gf2[:,end]

		return (gdl,gdb,gdr,gdt)
	elseif M.dim==3
		xf1,xf2,xf3 = getFaceGrids(M)
		
		gf1  = reshape(gd(xf1),tuple(M.n+[1,0,0]...))
		gx1   = gf1[1,:,:]
		gx2   = gf1[end,:,:]

		gf2  = reshape(gd(xf2),tuple(M.n+[0,1,0]...))
		gy1   = gf2[:,1,:]
		gy2   = gf2[:,end,:]

		gf3  = reshape(gd(xf3),tuple(M.n+[0,0,1]...))
		gz1   = gf3[:,:,1]
		gz2   = gf3[:,:,end]
		
		return (gx1,gy1,gz1,gx2,gy2,gz2)
	end
end