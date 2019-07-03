clc; clear;
domain = [0 2 0 1 0 1];
m      = [128, 64, 64];

W = zeros(m);

xc = reshape(getCellCenteredGrid(domain,m),[],3);
W = reshape((1.2*(xc(:,1)-.3).^2 + 2*(xc(:,2)-.6).^2 + .7*(xc(:,3)-.5).^2 ) < .1,m);
W = W + reshape((max(abs(xc-[.8 .2 .2]),[],2) < 0.1),m);
W = W + reshape((sum(abs(xc-[1 .7 .8]),2) < 0.1),m);

figure(1); clf;
volView(W,domain,m);
axis equal
axis(domain)
set(gca,'FontSize',20)
figDir = '/Users/lruthot/Dropbox/Projects/2018-MIPDECOpaper/images/3Dinstance';
 printFigure(gcf,fullfile(figDir,'3D-W.png'),'printFormat','.png','printOpts','-dpng');

save('Sources3D.mat','W','domain','m')

%%
load Sources3D.mat
%%
figure(2); clf;
figDir = '/Users/lruthot/Dropbox/Projects/2018-MIPDECOpaper/images/3Dinstance';
volView(utrue,domain,m+1,'isovalue',.06);
axis equal
hold on
for k=1:size(rec,1);
    plot3([rec(k,1); rec(k,1)],[rec(k,2); rec(k,2)],[0*rec(k,2); 0*rec(k,2)+1],'-r','LineWidth',2)
end
axis(domain)
set(gca,'FontSize',20)
printFigure(gcf,fullfile(figDir,'3D-utrue.png'),'printFormat','.png','printOpts','-dpng');

