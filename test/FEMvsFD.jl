
using PyPlot
using ConvDiff
using jInv.Mesh
using jInvVis
using jInv.LinearSolvers

domain      = [0 2 0 1]
nc		    = [128 64]*4
pNoise      = .05
noBoreHoles = 30
Mc          = getRegularMesh(domain,nc)
P,mf        = getRandomBoreholes2D(noBoreHoles,nc,nc) 		# Randomly selecting borehols
pForFEM     = getConvDiffFEMParam(Mc,sig=0.01,Ainv=getMUMPSsolver());
v = kron([1;0.0],ones(prod(Mc.n+2)))
pFor  = getConvDiffParam(Mc,v,sig=pForFEM.sig,Ainv=getMUMPSsolver(),bc=(:dir,:neu,:neu,:neu));
w = zeros(Mc.n[1],Mc.n[2]);
w[30:40,30:40] = 1;
viewImage2D(w,Mc)
title("source term")
colorbar()
pFor.sig

dobsFEM, = getData(vec(w),pForFEM);
dobs,    = getData(vec(w),pFor);

subplot(3,1,1)
viewImage2D(pForFEM.Fields,getPaddedMesh(Mc))
ylabel("FEM")
colorbar()
title("Fields")


subplot(3,1,2)
viewImage2D(pFor.Fields,Mc)
ylabel("FD")
colorbar()

subplot(3,1,3)
viewImage2D(getNodalAverageMatrix(Mc)*pForFEM.Fields-pFor.Fields,Mc)
colorbar()
ylabel("Difference")

println("relative error (PDE solve FEM vs FD): $(norm(vec(getNodalAverageMatrix(Mc)*pForFEM.Fields-pFor.Fields))/norm(vec(pFor.Fields)))")

"""
fdAnisoTV(W,Mc)

prototype implementation of regularizer using for loops 
"""
function fdAnisoTV(W,Mc)
    W = reshape(W,Mc.n[1], Mc.n[2])
    L = Mc.h
    
    Reg = 0.0;
    # compute first term of regularizer (derivatives in first dimension)
    for i=2:Mc.n[1]
        for j=1:Mc.n[2]
            Reg += Mc.h[2]*abs(W[i,j] - W[i-1,j])                
        end
    end
    # compute first term of regularizer (derivatives in second dimension)
    for i=1:Mc.n[1]
        for j=2:Mc.n[2]
            Reg += Mc.h[1]*abs(W[i,j] - W[i,j-1])                
        end
    end
    return Reg
end

Rfd = fdAnisoTV(w,Mc)
   
using jInv.InverseSolve
Rjinv = anisoTVReg(vec(w),0*vec(w),Mc,eps=0)[1]
println("Relative error (regularizer for loop): $(abs(Rfd-Rjinv)/abs(Rjinv))")

"""
fdPDE(U,Mc)

prototype implementation of PDE using for loops

here, assume that v = [1;0]
"""
function fdPDE(U,Mc;c=0.01)
    Up = zeros(Mc.n[1]+2, Mc.n[2]+2)
    Up[2:end-1,2:end-1] = reshape(U,Mc.n[1], Mc.n[2])
    
    # left boundary: 0-Dirichlet
    Up[1,2:end-1] = -U[1,:]
    # right boundary: 0-Neumann    
    Up[end,2:end-1] = U[end,:]
    # top boundary: 0-Neumann    
    Up[2:end-1,end] = U[:,end]
    # bottom boundary: 0-Neumann    
    Up[2:end-1,1] = U[:,1]

    L = Mc.h

    # allocate output (with ghostpoints for simpler indexing)
    V = zeros(Mc.n[1]+2, Mc.n[2]+2)
    for i=2:Mc.n[1]+1
        for j=2:Mc.n[2]+1
        V[i,j] = (c/(L[1]^2))*(2*Up[i,j]-Up[i-1,j]-Up[i+1,j]) + 
                 (c/(L[2]^2))*(2*Up[i,j]-Up[i,j-1]-Up[i,j+1]) +
                 (Up[i,j]-Up[i-1,j])/L[1]
        end
    end
    V = V[2:Mc.n[1]+1,2:Mc.n[2]+1]
    return V
end
U = randn(Mc.n[1],Mc.n[2])
# U = copy(w)
Vfd = fdPDE(copy(U),Mc,c=pFor.sig);

Vjinv = reshape(pFor.A*vec(U),Mc.n[1],Mc.n[2]);

subplot(3,1,1);
viewImage2D(Vfd,Mc)
colorbar()

subplot(3,1,2);
viewImage2D(Vjinv,Mc)
colorbar()

subplot(3,1,3);
viewImage2D(abs.(Vfd-Vjinv),Mc)
colorbar()

println("relative error (conv diff for loop): $(norm(vec(Vjinv-Vfd))/norm(vec(Vjinv)))")

