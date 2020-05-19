#####################################################
# PlotUW.ampl : plot the solutions (w and u) 
#####################################################
let leastSqr := (1/2/sigma)*sum{k in 1..K} (V[k]-b[k])^2;
let regularn := Lx*Ly*(sum{i in 2..Nx, j in 2..Ny} sqrt(gamma + 
((W[i,j]-W[i-1,j])/Lx)^2 + ((W[i,j]-W[i,j-1])/Ly)^2));
printf "\nAMPL output";
printf "\nalpha = %.5f, misfit term = %.5f, regularizer = %.5f", alpha, leastSqr, regularn;
##############################################################################
##### Generate Matlab code for printing ampl solution
##############################################################################
printf "\n\n%%% >>> Paste commands below in Matlab or Octave";
printf "\n%%% ========================================";
printf "\n%%%load instance file containing rec";
printf "\nNx = %f;\n", Nx;
printf "Ny = %f;\n", Ny;
printf "x1 = linspace(%f, %f, %d);\n", domain[1], domain[2], Nx;
printf "y1 = linspace(%f, %f, %d);\n", domain[3], domain[4], Ny;
printf "\n";

printf "W = [ "  ;
for {j in 1..Ny} {
printf{i in 1..Nx} "  %.5f\n", W[i,j]  ;
};
printf "  ];\n"  ;

printf "nfig = 1;\n";
printf "fig = figure(nfig); clf;\n"  ;
printf "nfig = nfig +1;\n"; 
printf "imagesc(x1, y1, (reshape(W,Nx,Ny))');\n";
printf "caxis([0 1]); hold on;\n";
printf "axis equal tight;\n";
printf "colorbar;\n";
printf "set(gca,'xtick',[]);\n";
printf "set(gca,'ytick',[]);\n";
printf "orient(fig,'landscape');\n"  ; 
printf "scatter(rec(:,1), rec(:,2), 30, 'filled', 'r');\n";
printf "ax = gca;\n";
printf "title('Sources (w) on %dx%d mesh');\n",Nx,Ny;
printf "ax.YDir = 'normal';\n";
printf "h=gcf;\n";
printf "set(h,'PaperOrientation','landscape');\n";
printf "set(h,'PaperUnits','normalized');\n";
printf "set(h,'PaperPosition', [0 0 1 1]);\n";
printf "print(gcf,'-dpdf','W%dx%d.pdf');\n\n\n",Nx,Ny;

#printf "U = [ "  ;
#for {j in 1..Ny+1} {
#printf{i in 1..Nx+1} "  %.5f\n", U[i,j]  ;
#}; 
#printf "  ];\n"  ;
#printf "fig = figure(nfig); clf;\n"  ;
#printf "nfig = nfig +1;\n"; 
#printf "imagesc(x1, y1, (reshape(U,(Nx+1),(Ny+1))'));\n";
#printf "caxis([0 1]); hold on;\n";
#printf "axis equal tight;\n";
#printf "colorbar;\n";
#printf "set(gca,'xtick',[]);\n";
#printf "set(gca,'ytick',[]);\n";
#printf "orient(fig,'landscape');\n"  ; 
#printf "scatter(rec(:,1), rec(:,2), 30, 'filled', 'r');\n";
#printf "ax = gca;\n";
#printf "ax.YDir = 'normal';\n";
#printf "h=gcf;\n";
#printf "set(h,'PaperOrientation','landscape');\n";
#printf "set(h,'PaperUnits','normalized');\n";
#printf "set(h,'PaperPosition', [0 0 1 1]);\n";
#printf "title('States (u) on %dx%d mesh');\n",Nx,Ny;
#printf "print(gcf,'-dpdf','U%dx%d.pdf');\n",Nx,Ny;
