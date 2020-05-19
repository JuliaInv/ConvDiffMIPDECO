######################################################################################
# SrcInvertFDM.mod - source inversion problem using finite-difference discretization
######################################################################################

######################################################################################
# Declare Parameters
param c default .01;                		# ... viscosity
param gamma default .001;           		# ... in regularization term
param alpha default 0.008531;       		# ... regularization parameter
param sigma default 1.0;            		# ... sigma estimate of variance
param domain{1..4};				# ... computational domain = [domain[1],domain[2]] x [domain[3],domain[4]]
param Nx >= 0, integer, default 4;  		# ... number of elements along x-axis
param Ny >= 0, integer, default 4;  		# ... number of elements along y-axis
param Lx:=(domain[2] - domain[1])/Nx;		# ... discretization step length along x-axis
param Ly:=(domain[4] - domain[3])/Ny;		# ... discretization step length along y-axis
param K > 0, integer;               		# ... number of measurements
param X{1..K};                      		# ... x coordinate of measurement locations
param Y{1..K};					# ... y coordinate of  measurement locations
param b{1..K};					# ... measurement values
param I{k in 1..K}  = floor(X[k]/Lx+0.5);       # ... derived expressions from (X,Y,b)
param J{k in 1..K}  = floor(Y[k]/Ly+0.5);   
param Xh{k in 1..K} = X[k] - I[k]*Lx + Lx/2;
param Yh{k in 1..K} = Y[k] - J[k]*Ly + Ly/2;
param leastSqr;					# ... misfit term in objective
param regularn;					# ... regularization term in objective

# Declare Variables
var U{0..Nx+1, 0..Ny+1};         		# ... recovered states
var W{1..Nx, 1..Ny} binary;         		# ... source locations
var V{k in 1..K}              	      		# ... interpolated U-value at (X[k],Y[k])
= (1/(Lx*Ly))*( U[I[k],J[k]]*(Lx-Xh[k])*(Ly-Yh[k]) + U[I[k],J[k]+1]*(Lx-Xh[k])*Yh[k] + 
U[I[k]+1,J[k]]*Xh[k]*(Ly-Yh[k]) + U[I[k]+1,J[k]+1]*Xh[k]*Yh[k]);


# Declare Objective 
minimize RegularizedLeastSqr: (1/2/sigma)*sum{k in 1..K} (V[k]-b[k])^2 +
alpha*Lx*Ly*(sum{i in 2..Nx, j in 2..Ny} sqrt(gamma + ((W[i,j]-W[i-1,j])/Lx)^2
+ ((W[i,j]-W[i,j-1])/Ly)^2));

# Declare Constraints 
subject to 

  # ... finite-difference approximation of advection-diffusion equation
  finiteDiff {i in 1..Nx, j in 1..Ny}: (c/(Lx*Ly))*(4*U[i,j]-U[i-1,j]-U[i+1,j]-U[i,j-1]-U[i,j+1]) + (U[i,j]-U[i-1,j])/Lx = W[i,j];

  # ... Dirichtlet boundary conditions on inflow boundary x=0
  boundaryConX0 {j in 0..Ny}: U[0,j] = -U[1,j]; 

  # ... Neuman boundary conditions
  boundaryConX1 {j in 0..Ny}: U[Nx+1,j] = U[Nx,j]; 
  boundaryConY0 {i in 0..Nx}: U[i,0] = U[i,1]; 
  boundaryConY1 {i in 0..Nx}: U[i,Ny+1] = U[i,Ny];
