export getFEMMatrices2D, getFEMMatrices3D, getFEMsource3D, getFEMsource2D, getFEMrhs3D, getFEMrhs2D

function Phi(i::Int,x::Float64,y::Float64)
    if i==1
        return (1-x)*(1-y)
    elseif i==2
        return  x*(1-y)
    elseif i==3
        return  x*y
    elseif i==4
        return y*(1-x)
    end 
end

function Phi(i::Int, x::Float64, y::Float64, z::Float64)
    if i==1
        return   (1-x)*(1-y)*(1-z)
    elseif i==2 
        return   x*(1-y)*(1-z)
    elseif i==3 
        return   x*y*(1-z)
    elseif i==4 
        return   y*(1-x)*(1-z)
	elseif i==5 
	    return   (1-x)*(1-y)*z
	elseif i==6 
	    return   x*(1-y)*z
	elseif i==7 
	    return   x*y*z
	elseif i==8 
	    return   y*(1-x)*z
    end 
end
function gradPhi(i::Int, x::Float64, y::Float64, z::Float64)
    if i==1
        return [-(1-y)*(1-z) -(1-x)*(1-z) -(1-x)*(1-y)]
    elseif i==2
        return [(1-y)*(1-z) -x*(1-z) -x*(1-y)]
    elseif i==3
        return [y*(1-z) x*(1-z) -x*y]
    elseif i==4
        return [-y*(1-z) (1-x)*(1-z) -y*(1-x)]
    elseif i==5
        return [-(1-y)*z -(1-x)*z  (1-x)*(1-y)]
    elseif i==6
        return [(1-y)*z -x*z  x*(1-y)]
    elseif i==7
        return [y*z x*z x*y]
    elseif i==8
        return [-y*z (1-x)*z y*(1-x)]
    end 
end

function gradPhi(i::Int,x::Float64, y::Float64)
    if i==1
        return  [y-1 x-1]
    elseif i==2
        return  [1-y -x]
    elseif i==3
        return  [y x]
    elseif i==4
        return  [-y 1-x]
    end 
end




function getFEMMatrices2D(M::RegularMesh;sig=0.01,v=[1.;0])
h      = M.h;
Jac    = Diagonal(h)
invJac = inv(Jac)
detJac = prod(h)

gaussRule=5
wGauss = [0.118463442528094; 0.239314335249683; 0.284444444444445;
0.239314335249683; 0.118463442528094];
xGauss = [0.9530899229693319; 0.7692346550528416; 0.5000000000000000;
0.2307653449471585; 0.0469100770306681];
xGauss = sort(xGauss)


X,Y = meshgrid(xGauss,xGauss);
W   = kron(wGauss,wGauss);
# wGauss*wGauss';
numquadnodes = gaussRule^M.dim;
X = vec(X)
Y = vec(Y)
W = vec(W)

# element stiffness matrix for Laplacian
sm_laplace = zeros(4,4)
for ii=1:4
    for jj=1:4
        for q=1:numquadnodes
            x = X[q];
            y = Y[q];
            sm_laplace[ii,jj] +=  dot(gradPhi(ii,x,y)*invJac,gradPhi(jj,x,y)*invJac)*W[q]*detJac;
        end
    end
end

# element mass matrix
mass = zeros(4,4);
for ii=1:4
    for jj=1:4
        for q=1:numquadnodes
            x = X[q];
            y = Y[q];
            mass[ii,jj] += Phi(ii,x,y)*Phi(jj,x,y)*W[q]*detJac
        end
    end
end

# element mass matrix with constant ansatz functions
mass_const = zeros(4);
for ii=1:4
    for q=1:numquadnodes
        x = X[q];
        y = Y[q];
        mass_const[ii] +=  Phi(ii,x,y)*1.0*W[q]*detJac
    end
end

# element stiffness matrix for advection
sm_advection = zeros(4,4)
for ii=1:4
    for jj=1:4
        for q=1:numquadnodes
            x = X[q]
            y = Y[q]
            sm_advection[ii,jj] +=  Phi(ii,x,y)*dot(v,gradPhi(jj,x,y)*invJac)*W[q]*detJac
        end
    end
end


nn = M.n+1;

Ism = zeros(Int,0)
Jsm = zeros(Int,0)
Vsm = zeros(Float64,0)

Imass = zeros(Int,0)
Jmass = zeros(Int,0)
Vmass = zeros(Float64,0)

Imass_const = zeros(Int,0)
Jmass_const = zeros(Int,0)
Vmass_const = zeros(Float64,0)

for ii=1:nn[2]-1
    for jj=1:nn[1]-1
        nodes = [jj+nn[1]*(ii-1); jj+nn[1]*(ii-1)+1; jj+nn[1]*ii+1; jj+nn[1]*ii];
        for k=1:4
            for l=1:4
                push!(Imass,nodes[k])
				push!(Jmass,nodes[l])
				push!(Vmass,mass[k,l])
				
				push!(Ism, nodes[k])
				push!(Jsm, nodes[l])
				vsm = sig*sm_laplace[k,l] + sm_advection[k,l]
				push!(Vsm, vsm)
           end
			push!(Imass_const,nodes[k])
			push!(Jmass_const,(ii-1)*(nn[1]-1)+jj)
			push!(Vmass_const,mass_const[k])
        end
    end
end
SM         = sparse(Ism,Jsm,Vsm, prod(nn), prod(nn))
Mass       = sparse(Imass, Jmass, Vmass, prod(nn), prod(nn))
Mass_const = sparse(Imass_const, Jmass_const, Vmass_const, prod(nn), M.nc)
return Mass, Mass_const, SM
end

function getFEMMatrices3D(M::RegularMesh;sig=0.01,v=[1.;0;0])

h      = M.h;
Jac    = Diagonal(h)
invJac = inv(Jac)
detJac = prod(h)

gaussRule=5
wGauss = [0.118463442528094; 0.239314335249683; 0.284444444444445;
0.239314335249683; 0.118463442528094];
xGauss = [0.9530899229693319; 0.7692346550528416; 0.5000000000000000;
0.2307653449471585; 0.0469100770306681];

X,Y,Z = meshgrid(xGauss,xGauss,xGauss);
W     = kron(wGauss,kron(wGauss,wGauss));
# wGauss*wGauss';
numquadnodes = gaussRule^M.dim;
X = vec(X)
Y = vec(Y)
Z = vec(Z)
W = vec(W)

# element stiffness matrix for Laplacian
sm_laplace = zeros(8,8)
for ii=1:8
    for jj=1:8
        for q=1:numquadnodes
            x = X[q]
            y = Y[q]
			z = Z[q]
			sm_laplace[ii,jj] +=  dot(gradPhi(ii,x,y,z)*invJac,gradPhi(jj,x,y,z)*invJac)*W[q]*detJac;
        end
    end
end

# element mass matrix
mass = zeros(8,8);
for ii=1:8
    for jj=1:8
        for q=1:numquadnodes
            x = X[q]
            y = Y[q]
            z = Z[q]
			mass[ii,jj] += Phi(ii,x,y,z)*Phi(jj,x,y,z)*W[q]*detJac
        end
    end
end

# element mass matrix with constant ansatz functions
mass_const = zeros(8);
for ii=1:8
    for q=1:numquadnodes
        x = X[q]
        y = Y[q]
        z = Z[q]
		mass_const[ii] +=  Phi(ii,x,y,z)*1.0*W[q]*detJac
    end
end

# element stiffness matrix for advection for v = [1;0;0]
sm_advection = zeros(8,8)
for ii=1:8
    for jj=1:8
        for q=1:numquadnodes
            x = X[q]
            y = Y[q]
            z = Z[q]
			sm_advection[ii,jj] += Phi(ii,x,y,z)*dot(v,gradPhi(jj,x,y,z)*invJac)*W[q]*detJac
        end
    end
end


nn         = M.n+1;


Ism = zeros(Int,0)
Jsm = zeros(Int,0)
Vsm = zeros(Float64,0)

Imass = zeros(Int,0)
Jmass = zeros(Int,0)
Vmass = zeros(Float64,0)

Imass_const = zeros(Int,0)
Jmass_const = zeros(Int,0)
Vmass_const = zeros(Float64,0)

nodes = zeros(Int,8)
incz = prod(nn[1:2])

for kk=1:nn[3]-1
	for ii=1:nn[2]-1
	    for jj=1:nn[1]-1
	 		nodes[1] =  jj+nn[1]*(ii-1)+incz*(kk-1);
	 		nodes[2] =  jj+nn[1]*(ii-1)+1+incz*(kk-1); 
			nodes[3] =  jj+nn[1]*ii+1+incz*(kk-1); 
			nodes[4] =  jj+nn[1]*ii+incz*(kk-1)
			nodes[5] =  jj+nn[1]*(ii-1)+incz*kk; 
			nodes[6] =  jj+nn[1]*(ii-1)+1+incz*kk; 
			nodes[7] =  jj+nn[1]*ii+1+incz*kk; 
			nodes[8] =  jj+nn[1]*ii+incz*kk;
	       for k=1:8
	            for l=1:8
	                push!(Imass,nodes[k])
					push!(Jmass,nodes[l])
					push!(Vmass,mass[k,l])
					
					push!(Ism, nodes[k])
					push!(Jsm, nodes[l])
					vsm = sig*sm_laplace[k,l] + sm_advection[k,l]
					push!(Vsm, vsm)
	           end
			end
			for k=1:8
				push!(Imass_const,nodes[k])
				push!(Jmass_const,(ii-1)*(nn[1]-1)+jj+ prod(nn[1:2]-1)*(kk-1))
				push!(Vmass_const,mass_const[k])
	        end
	    end
	end
end
SM         = sparse(Ism,Jsm,Vsm, prod(nn), prod(nn))
Mass       = sparse(Imass, Jmass, Vmass, prod(nn), prod(nn))
Mass_const = sparse(Imass_const, Jmass_const, Vmass_const, prod(nn), M.nc)
return Mass, Mass_const, SM
end


function getFEMsource3D(M::RegularMesh)
	xc  = getCellCenteredGrid(M)
	blob1 = (xc[:,1].>.1) .& (xc[:,1] .< 2.2 ) .& (xc[:,2].>.3) .& (xc[:,2] .< 0.5) .& (xc[:,3] .> .3 ) .& (xc[:,3] .< 1.1);
	f     = 1.0*(blob1 )
	return f
end

function getFEMsource2D(M::RegularMesh)

	xc  = getCellCenteredGrid(M)
	blob1 = (xc[:,1].>.5) .& (xc[:,1] .< .7 ) .& (xc[:,2].>.3) .& (xc[:,2] .< 0.5)
	blob2 = (xc[:,1].>1.0) .& (xc[:,1] .< 1.5 ) .& (xc[:,2].>.6) .& (xc[:,2] .< 0.8)
	f     = 1.0*(blob1 + blob2)
	return f
end


function getFEMrhs2D(M::RegularMesh)
nn  = M.n 
hh  = M.h 
rhs = zeros(M.nc,1);
for ii=1:nn[2]
	for jj=1:nn[1]
		x      = (jj-0.5)*hh[1];
		y      = (ii-0.5)*hh[2];
		blob1  = (x > 0.5) & (x < 0.7 ) & (y > 0.3) & (y < 0.5)
		blob2  = (x > 1.0) & (x < 1.5 ) & (y > 0.6) & (y < 0.8)
		p      =  (ii-1)*nn[1]+jj
		rhs[p] = 1.0*(blob1 + blob2)
		
		
	end 
end
return rhs 
end


function getFEMrhs3D(M::RegularMesh)
nn   = M.n 
hh   = M.h 
rhs  = zeros(M.nc,1)
incz = prod(nn[1:2])
for kk=1:nn[3]
	for ii=1:nn[2]
		for jj=1:nn[1]
			x      = (jj-0.5)*hh[1];
			y      = (ii-0.5)*hh[2];
			z      = (kk-0.5)*hh[3];
			blob1  = (x >.1) & (x < 1 ) & (y > 0.3) & (y < 0.5) & (z > 0.3 ) & (z < .5);
			blob2  = (x > 1) & (x < 2 ) & (y > 0.6) & (y < 0.9) & (z > 0.6 ) & (z < .9);
			p      = (ii-1)*nn[1]+jj+(kk-1)*incz 
			#rhs[p] = 1.0*(blob1)
			rhs[p] = 1.0*(blob1 + blob2)
		end 
	end 
end
return rhs 
end

# 
# xnodes = linspace(0,xlength, nn[1]);
# ynodes = linspace(0,ylength, nn[2]);
# [X,Y] = meshgrid(xnodes, ynodes);
# 
# rhs = zeros((nn[1]-1)*(nn[2]-1),1);
# for ii=1:(nn[2]-1)
#     for jj=1:(nn[1]-1)
#         x = (jj+0.5)*hx;
#         y = (ii+0.5)*hy;
#         rhs((ii-1)*(nn[1]-1)+jj) = (x>0.5 && x<0.7 && y>0.3 && y<0.5) || (x>1.0 && x<1.5 && y>0.6 && y<0.8);
#     end
# end
# V = Mass_const*rhs;
# u = SM\V;                         % ... solve SM*u = Mass_const*rhs, where whs are element controls
#                                   % ... u has Nx*nn[2] entries U(i,j) = u((i-1)*Nx+j) for i=1,Nx, j=1,nn[2]
#                                   % ... rhs has (Nx-1)*(Ny-1) entries W(i,j) = rhs((i-1)*(Nx-1)+j) for i=1,Nx-1, j=1,Ny-1
# PrintAMPL = 0;
# if PrintAMPL
#   outFid = fopen('FEM.dat','w');
#   for k=1:Nx*Ny
#     % SM(k,:)*u = Mass_const(k,:)*rhs
#     irows = find(SM(k,:));
#     fprintf(outFid,'FEM_%-6i:',k);
#     for l=1:length(irows)
#       jj = mod(irows(l),Nx) + 1;
#       ii = floor(irows(l) / Nx) + 1;
#       aa = full(SM(k,irows(l)));
#       fprintf(outFid,' + %12.8g * U[%6i,%6i]', aa,ii,jj);
#     end
#     fprintf(outFid,';\n');
#   end
#   fclose(outFid);
# end
# 
# u_plot = reshape(u, Nx, Ny)';
# figure(1); clf;
# imagesc(reshape(rhs,Nx-1,Ny-1))
# figure(2); clf;
# surf(X,Y,u_plot)
# %daspect("auto")
# pbaspect([xlength ylength max(u)-min(u)])
