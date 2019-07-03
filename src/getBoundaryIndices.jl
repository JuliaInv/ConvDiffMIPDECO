export getBoundaryIndices, getBoundaryIndicesFEM


"""
function getBoundaryIndices(M::RegularMesh, bc=(:dir,:neu,:neu,:neu))

gets indices of boundary points, neighbors of boundary points and interior
points in extended mesh

Input: 

	M::RegularMesh - mesh
	bc             - chooses Dirichlet or Neuman BC for each boundary
	                 Ex1: dim=2, domain=[0,1,0,1] --> bc is a vector of length 4
					      edge 1: x1=0, x2 in [0,1] (left)
						  edge 2: x2=0, x1 in [0,1] (bottom)
						  edge 3: x1=1, x2 in [0,1] (right)
						  edge 4: x2=1, x1 in [0,1] (top)
		             Ex1: dim=3, domain=[0,1,0,1,0,1] --> bc is a vector of length 6
		 			      face 1: x1=0, x2,x3 in [0,1]
		 			 	  face 2: x2=0, x1,x3 in [0,1] 
		 			 	  face 3: x3=0, x1,x2 in [0,1] 
		 			      face 4: x1=1, x2,x3 in [0,1]
		 			 	  face 5: x2=2, x1,x3 in [0,1] 
		 			 	  face 6: x3=3, x1,x2 in [0,1] 

Outputs:
	iddir     - indices of ghost points associated with Dirichlet conditions				  
	idneu     - indices of ghost points associated with Neuman conditions				  
	iddirn    - indices of neighbors of ghost points associated with Dirichlet conditions				  
	idneun    - indices of neighbors of ghost points associated with Neuman conditions				  
	idint     - indices of interior nodes	
				   
"""
function getBoundaryIndices(M::RegularMesh, bc=(:dir,:neu,:neu,:neu,:neu,:neu))
	
	# get boundary conditions
	idg,idgn,idint = getGhostIndices(M.n)
	iddir  = zeros(Int,0)
	idneu  = zeros(Int,0)
	iddirn = zeros(Int,0)
	idneun = zeros(Int,0)
	for k=1:length(bc)
		if bc[k]==:dir
			iddir  = vcat(iddir,idg[k])
			iddirn = vcat(iddirn, idgn[k])
		elseif bc[k]==:neu
			idneu  = vcat(idneu,idg[k])
			idneun = vcat(idneun,idgn[k])
		end
	end
	
	return iddir, idneu, iddirn, idneun, idint
end


"""
function getBoundaryIndicesFEM(M::RegularMesh, bc=(:dir,:neu,:neu,:neu))

Input: 

	M::RegularMesh - mesh
	bc             - chooses Dirichlet or Neuman BC for each boundary
	                 Ex1: dim=2, domain=[0,1,0,1] --> bc is a vector of length 4
					      edge 1: x1=0, x2 in [0,1] (left)
						  edge 2: x2=0, x1 in [0,1] (bottom)
						  edge 3: x1=1, x2 in [0,1] (right)
						  edge 4: x2=1, x1 in [0,1] (top)
		             Ex1: dim=3, domain=[0,1,0,1,0,1] --> bc is a vector of length 6
		 			      face 1: x1=0, x2,x3 in [0,1]
		 			 	  face 2: x2=0, x1,x3 in [0,1] 
		 			 	  face 3: x3=0, x1,x2 in [0,1] 
		 			      face 4: x1=1, x2,x3 in [0,1]
		 			 	  face 5: x2=2, x1,x3 in [0,1] 
		 			 	  face 6: x3=3, x1,x2 in [0,1] 

Outputs:
	iddir     - indices of FEM nodes associated with Dirichlet conditions				  
				   
"""
function getBoundaryIndicesFEM(M::RegularMesh, bc=(:dir,:neu,:neu,:neu,:neu,:neu))
	
	# get all nodes
	ids    = reshape(1:prod(M.n.+1),tuple(M.n.+1...))
	iddir  = zeros(Int,0)
	
	if M.dim==2
		if bc[1]==:dir 
			iddir = vcat(iddir,vec(ids[1,:]));
		end
		if bc[2]==:dir 
			iddir = vcat(iddir,vec(ids[:,1]));
		end
		if bc[3]==:dir 
			iddir = vcat(iddir,vec(ids[end,:]));
		end
		if bc[4]==:dir 
			iddir = vcat(iddir,vec(ids[:,end]));
		end
	elseif M.dim==3
		if bc[1]==:dir 
			iddir = vcat(iddir,vec(ids[1,:,:]));
		end
		if bc[2]==:dir 
			iddir = vcat(iddir,vec(ids[:,1,:]));
		end
		if bc[3]==:dir 
			iddir = vcat(iddir,vec(ids[:,:,1]));
		end
		if bc[4]==:dir 
			iddir = vcat(iddir,vec(ids[end,:,:]));
		end
		if bc[5]==:dir
			iddir = vcat(iddir,vec(ids[:,end,:]));
		end
		if bc[6]==:dir 
			iddir = vcat(iddir,vec(ids[:,:,end]));
		end
	end
	return unique(iddir)
end
