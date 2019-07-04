export getLaplacianMatrix, getConvectionMatrix

"""
function ddxc(n)

centered one-dimensional finite difference operator.

Input: 

    n  - number of cells

Output:

    dx - sparse matrix for 1D derivative size(dx) = (n,n)
"""
function ddxc(n::Int)
    dx = spdiagm(-1=>fill(-1.0,n-1), 1=>fill(1.0,n-1))    
    return dx
end

"""
function ddxfwd(n)

forward one-dimensional finite difference operator.

Input: 

    n  - number of cells

Output:

    dx - sparse matrix for 1D derivative size(dx) = (n,n)
"""
function ddxfwd(n::Int)
    dx = spdiagm(0=>fill(-1.0,n),1=>fill(1.0,n-1))    
    return dx
end

"""
function ddxbwd(n)

backward one-dimensional finite difference operator.

Input: 

    n  - number of cells

Output:

    dx - sparse matrix for 1D derivative size(dx) = (n,n)
"""
function ddxbwd(n::Int)
    dx = spdiagm(-1=>fill(-1.0,n-1),0=>fill(1.0,n))    
    return dx
end

"""
function getLaplacian(M::RegularMesh)

builds 2D/3D Laplacian for padded cell-centered vector including
ghost points and corners.

Input:

	M  - computational mesh
	
Output:

	Lap - 2D/3D Laplacian

"""
function getLaplacianMatrix(M::RegularMesh)
    if M.dim==2
		dx = Mesh.ddx(M.n[1]+1)/M.h[1]
		dy = Mesh.ddx(M.n[2]+1)/M.h[2]
		
		Dxx = kron(sparse(1.0I,M.n[2]+2, M.n[2]+2),dx'*dx)
		Dyy = kron(dy'*dy, sparse(I, M.n[1]+2, M.n[1]+2))
		return Dxx+Dyy
	elseif M.dim==3
		dx = Mesh.ddx(M.n[1]+1)/M.h[1]
		dy = Mesh.ddx(M.n[2]+1)/M.h[2]
		dz = Mesh.ddx(M.n[3]+1)/M.h[3]
		
		Dxx = kron(sparse(1.0I, M.n[3]+2, M.n[3]+2),kron(sparse(1.0I, M.n[2]+2, M.n[2]+2),dx'*dx))
		Dyy = kron(sparse(1.0I, M.n[3]+2, M.n[3]+2),kron(dy'*dy, sparse(1.0I, M.n[1]+2, M.n[1]+2)))
		Dzz = kron(dz'*dz, kron(sparse(1.0I, M.n[2]+2, M.n[2]+2), sparse(1.0I, M.n[1]+2, M.n[1]+2)))
		
		return Dxx + Dyy + Dzz
	end
end

"""
function getConvectionMatrix(M,v)
	
builds stationary convection operator on extended mesh

2D: C(v_1,v_2)     = Diagonal(v_1)*Dx + Diagonal(v_2)*Dy

3D: C(v_1,v_2,v_3) = Diagonal(v_1)*Dx + Diagonal(v_2)*Dy + Diagonal(v_3)*Dz



Input:

	M :: RegularMesh - original domain (not padded)
	v                - velocities discretized on cell-centers of extended mesh
	 				   or function that returns values (X -> vFun) where 
					   X has x,y(,z) coordinates column wise.
	
Output:

	Conv - 2D Convection operator

"""
function getConvectionMatrix(M::RegularMesh,v::Array{Float64})

	n   = M.n
	h   = M.h
	dim = M.dim
	v  = reshape(v,prod(n.+2),dim)

	if dim==2
    	dxf = ddxfwd(n[1]+2)./(h[1])
    	dxb = ddxbwd(n[1]+2)./(h[1])
    	dyf = ddxfwd(n[2]+2)./(h[2])
    	dyb = ddxbwd(n[2]+2)./(h[2])
    	
    	Dxf = kron(sparse(1.0I, n[2]+2, n[2]+2), dxf)
    	Dxb = kron(sparse(1.0I, n[2]+2, n[2]+2), dxb)
    	Dyf = kron(dyf, sparse(1.0I, n[1]+2, n[1]+2))
    	Dyb = kron(dyb, sparse(1.0I, n[1]+2, n[1]+2))
    	return Diagonal(max.(v[:,1],0))*Dxb + Diagonal(min.(v[:,1],0))*Dxf + 
		       Diagonal(max.(v[:,2],0))*Dyb + Diagonal(min.(v[:,2],0))*Dyf
	elseif dim==3
		dxf = ddxfwd(n[1]+2)./(M.h[1])
		dxb = ddxbwd(n[1]+2)./(M.h[1])
		dyf = ddxfwd(n[2]+2)./(M.h[2])
		dyb = ddxbwd(n[2]+2)./(M.h[2])
		dzf = ddxfwd(n[3]+2)./(M.h[3])
		dzb = ddxbwd(n[3]+2)./(M.h[3])
    
		Dxf = kron(sparse(1.0I, n[3]+2, n[3]+2),kron(sparse(1.0I, n[2]+2, n[2]+2), dxf));
		Dxb = kron(sparse(1.0I, n[3]+2, n[3]+2),kron(sparse(1.0I, n[2]+2, n[2]+2), dxb));
		Dyf = kron(sparse(1.0I, n[3]+2, n[3]+2),kron(dyf, sparse(1.0I, n[1]+2, n[1]+2)))
		Dyb = kron(sparse(1.0I, n[3]+2, n[3]+2),kron(dyb, sparse(1.0I, n[1]+2, n[1]+2)))
		Dzf = kron(dzf, sparse(1.0I, prod(n[1:2].+2), prod(n[1:2].+2)))
		Dzb = kron(dzb, sparse(1.0I, prod(n[1:2].+2), prod(n[1:2].+2)))
    	return Diagonal(max.(v[:,1],0))*Dxb + Diagonal(min.(v[:,1],0))*Dxf + 
               Diagonal(max.(v[:,2],0))*Dyb + Diagonal(min.(v[:,2],0))*Dyf +
               Diagonal(max.(v[:,3],0))*Dzb + Diagonal(min.(v[:,3],0))*Dzf 
	end
end
	
function getConvectionMatrix(M::RegularMesh,v::Function)
	xc = getNodalGrid(getPaddedMesh(M))
	return getConvectionMatrix(M,v(xc))
end