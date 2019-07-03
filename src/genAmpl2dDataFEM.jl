export genAmpl2dDataFEM
function
genAmpl2dDataFEM(M::RegularMesh,SM::SparseMatrixCSC{Float64},MC::SparseMatrixCSC{Float64},b::Vector,P::SparseMatrixCSC,dobs::Vector,nbx::Vector,nby::Vector,pNoise::Int,ns::Int)

#-------------------------------------------------------------------------------
# Generating lookup table for variables u and w
#-------------------------------------------------------------------------------
  nn1  = M.n+1 
  nn2 = nn1-1 
  np1 = prod(nn1)
  np2 = prod(nn2)
  wx  = zeros(Int64,np2,1)
  wy  = zeros(Int64,np2,1)
  ux  = zeros(Int64,np1,1)
  uy  = zeros(Int64,np1,1)
  #bc=(:dir,:neu,:neu,:neu,:neu,:neu)
  #iddir = getBoundaryIndicesFEM(M, bc)  # Dirichlet nodes
  #(dirX, dirY) = ind2sub(zeros(nn1[1],nn1[2]),iddir)
  for ii=1:nn2[2]
  	for jj=1:nn2[1]
  	wx[(ii-1)*nn2[1]+jj] = jj; 
  	wy[(ii-1)*nn2[1]+jj] = ii;
  	nodes =	[jj+(nn2[1]+1)*(ii-1);jj+(nn2[1]+1)*(ii-1)+1;jj+(nn2[1]+1)*(ii)+1;jj+(nn2[1]+1)*ii];
  	ux[nodes[1]] = jj;  uy[nodes[1]] = ii;
  	ux[nodes[2]] = jj+1; uy[nodes[2]] = ii;
  	ux[nodes[3]] = jj+1; uy[nodes[3]] = ii+1;
  	ux[nodes[4]] = jj;   uy[nodes[4]] = ii+1;		  
    end 
  end

#-------------------------------------------------------------------------------
# Generating data file for AMPL
#-------------------------------------------------------------------------------

	printAMPL = true
	fname = string(M.n[1], "x",	M.n[2],"S",ns,"P",length(dobs),"N",pNoise,".dat"); 	
	fname1 = string(M.n[1], "x",M.n[2],"S",ns,"P",length(dobs),"N",pNoise); 	
	if printAMPL
		open(fname,"w") do f
			println(f,"data;")
			println(f, "param outfile:=",fname1, ";")
			println(f,"param Ly  := $(M.domain[2]-M.domain[1]) ;")
			println(f,"param Lx  := $(M.domain[4]-M.domain[3]) ;")
			println(f,"param Ny  := $((M.n[1])+1) ;")
			println(f,"param Nx  := $((M.n[2])+1) ;")
			println(f,"param NBH := $(length(dobs)) ;")
			@printf(f,"param nbx := ")
			for k=1:length(dobs)
				@printf(f, " %d %.8f",k, nbx[k]) 
			end	
			@printf(f, ";\n") 
			@printf(f,"param nby := ")
			for k=1:length(dobs)
				@printf(f, " %d %.8f",k, nby[k]) 
			end	
			@printf(f, ";\n") 
			for k=1:length(dobs)
				irow = find(P[k,:])
				for i=1:length(irow)
					@printf(f, "let Uhat[%d,%d]:=%.8f;\n",ux[irow[i]],uy[irow[i]],dobs[k]) 
				end
			end	
			#for k=1:length(dirX)
			  #@printf(f, "fix U[%d,%d]:=0;\n",dirX[k], dirY[k]) 
			#end	
			for k=1:length(dobs)
				irow = find(P[k,:])
				for i=1:length(irow)
					@printf(f, "let P[%d,%d]:=%d;\n",ux[irow[i]],uy[irow[i]],1) 
				end
			end	

			println(f, "model;")
			println(f, "subject to")
			for k=1:nn1[1]*nn1[2]
			  #if mod(k,nn1[1]) != 1
			    irows = find(SM[k,:])
			    print(f, "FEM_",k,": ")
			    for i=1:length(irows)
				  #@printf(f,"+ U[%d,%d]*SM[%d,%d]",ux[irows[i]],uy[irows[i]],k,irows[i])
				  if i==1
					@printf(f,"%.8f* U[%d,%d]", SM[k,irows[i]],ux[irows[i]],uy[irows[i]])

				  else
					@printf(f,"+ %.8f* U[%d,%d]", SM[k,irows[i]],ux[irows[i]],uy[irows[i]])

				  end
			    end
			    print(f,"=")
			    irows = find(MC[k,:])
			    for i=1:length(irows)
				  #@printf(f,"+ W[%d,%d]*MC[%d,%d]",wx[irows[i]],wy[irows[i]],k,irows[i])
				  if i==1
					@printf(f,"%.8f* W[%d,%d] + %.8f", MC[k,irows[i]],wx[irows[i]],wy[irows[i]], b[k])
				  else
					@printf(f,"+ %.8f* W[%d,%d] + %.8f", MC[k,irows[i]],wx[irows[i]],wy[irows[i]], b[k])
				  end
			    end
				if length(irows) == 0 
			      @printf(f,"0")
				end
			    println(f,";")
			  #end
			end
		end
	end
	return
end

