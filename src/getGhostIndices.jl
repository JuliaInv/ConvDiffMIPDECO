export getGhostIndices

"""
function getGhotIndiced(n)
	
get indices of ghost points, boundary points and interior points for 
cell-centered grid with n cells

Input: 

   n::Array{Int}    -  number of cells
   
Output:

   idg   - indices of ghost points
   idgn  - indices of ghost point neighbors
   idint - indices of interior nodes
   
"""
function getGhostIndices(n::Array{Int})
	dim = length(n)
	
	if dim==2
		idg   = Array{Array{Int,1}}(4)
		idgn  = Array{Array{Int,1}}(4)
		
		idn    = reshape(collect(1:prod(n+2)),tuple(n+2...))::Array{Int,2}
		idg[1] = vec(idn[1,2:end-1])
		idg[2] = vec(idn[2:end-1,1])
		idg[3] = vec(idn[end,2:end-1])
		idg[4] = vec(idn[2:end-1,end])
		
		idgn[1] = vec(idn[2,2:end-1])
		idgn[2] = vec(idn[2:end-1,2])
		idgn[3] = vec(idn[end-1,2:end-1])
		idgn[4] = vec(idn[2:end-1,end-1])
		idint   = vec(idn[2:end-1,2:end-1])
	else
		idg   = Array{Array{Int,1}}(6)
		idgn  = Array{Array{Int,1}}(6)

		idn   = reshape(collect(1:prod(n+2)),tuple(n+2...))::Array{Int,3}
		idg[1] = vec(idn[1,2:end-1,2:end-1])
		idg[2] = vec(idn[2:end-1,1,2:end-1])
		idg[3] = vec(idn[2:end-1,2:end-1,1])
		idg[4] = vec(idn[end,2:end-1,2:end-1])
		idg[5] = vec(idn[2:end-1,end,2:end-1])
		idg[6] = vec(idn[2:end-1,2:end-1,end])
		
		idgn[1] = vec(idn[2,2:end-1,2:end-1])
		idgn[2] = vec(idn[2:end-1,2,2:end-1])
		idgn[3] = vec(idn[2:end-1,2:end-1,2])
		idgn[4] = vec(idn[end-1,2:end-1,2:end-1])
		idgn[5] = vec(idn[2:end-1,end-1,2:end-1])
		idgn[6] = vec(idn[2:end-1,2:end-1,end-1])
		
		idint = vec(idn[2:end-1,2:end-1,2:end-1])
	end
	return idg,idgn,idint
end
