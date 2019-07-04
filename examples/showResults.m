clear;

load('Peaks2D.mat')
dataset = 'Peaks2D';
reg     = 'wTVReg';
m       = [256 128];
domain = [0 2 0 1];
doPrint = 1;


files = dir([dataset '-' reg '-noise-0.1-' num2str(m(1)) 'x' num2str(m(2)) '.mat']);
figDir = '/Users/lruthot/Dropbox/Projects/2019-MIPDECOpaper/images/figPeaks2D/';



fname = files(1).name;
load(fname)
%%
idx = 11;
fig = figure();
fig.Name = fname;
loglog(Reg,Mis,'-x','LineWidth',2,'MarkerSize',8)
hold on
loglog(Reg([1;idx;30]),Mis([1;idx;30]),'or','MarkerSize',10)
ylabel('misfit');
xlabel('regularizer');
legend()
if doPrint
      matlab2tikz(fullfile(figDir,regexprep(fname,'.mat','.tex')),'width','\iwidth','height','\iheight');
end


%%
fig = figure();
fig.Name = fname;
src = flip(permute(reshape(Sources,m(1),m(2),[]),[2 1 3]),1);
montageArray(src);
axis equal tight

%%
src = reshape(Sources,m(1),m(2),[]);
for k=1:size(src,3);
    fig = figure(); clf;
    viewImage2Dsc(src(:,:,k),domain,m)
    caxis([0 1]);
    axis off
    if doPrint
        printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-' num2str(k) '.png']),'printOpts','-dpng','printFormat','.png');
        close(fig);
    end
end

%% plot solution of MIPDECO
xa = linspace(domain(1),domain(2),size(W,1));
ya = linspace(domain(3),domain(4),size(W,2));
[XA,YA] = meshgrid(xa,ya);
xfine =getCellCenteredGrid(domain,size(W));

for k=1:size(MIPDECO,3)
    for j=1:size(MIPDECO,4)
        fig = figure(); clf;
        viewImage2Dsc(MIPDECO(:,:,k,j),domain,m);
        hold on;
        contour(XA',YA',W,1,'-r','LineWidth',3);
        caxis([0 1]);
        axis off
        
        mip = reshape(nnInter(MIPDECO(:,:,k,j),domain,xfine),size(W));
        intersect = mip.*W;
        union     = (mip+W)>eps;
        IoU = nnz(intersect)/nnz(union);
        
        
        text(1.2,.1,sprintf('IoU=%1.1f%%',IoU*100),'FontSize',40,'Color',[.99 .99 .99])
        if doPrint
            
            printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-MIPDECO-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
            close(fig);
        end
    end
end

%% plot initiliztion of MIPDECO
xa = linspace(domain(1),domain(2),size(W,1));
ya = linspace(domain(3),domain(4),size(W,2));
[XA,YA] = meshgrid(xa,ya);

for k=1:size(MIPDECO,3)
    for j=1:size(MIPDECO,4)
        fig = figure(); clf;
        viewImage2Dsc(MIPDECOinit(:,:,k,j),domain,m);
        hold on;
        contour(XA',YA',W,1,'-r','LineWidth',3);
        caxis([0 1]);
        
        mip = reshape(nnInter(MIPDECOinit(:,:,k,j),domain,xfine),size(W));
        intersect = mip.*W;
        union     = (mip+W)>eps;
        IoU = nnz(intersect)/nnz(union);
        
        
        text(1.2,.1,sprintf('IoU=%1.1f%%',IoU*100),'FontSize',40,'Color',[.99 .99 .99])
        
        
        axis off
        if doPrint
            printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-MIPDECOinit-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
            close(fig);
        end
    end
end
%%


%%
if doPrint
   for k=1:size(His,3)
       for j=1:size(His,4)
           his = His(:,:,k,j)
           fig = figure(); clf;
           iter = numel(find(his(:,1)));
           yyaxis left
           plot(his(1:iter,1),'linewidth',3)
           ylabel('objective function')
           
           yyaxis right
           plot(his(1:iter,3),'linewidth',3)
           ylabel('Trust Region radius')
           xlabel('iterations')
           set(gca,'FontSize',20)
           axis([0 50 0 512]);
           printFigure(fig,fullfile(figDir,[dataset '-' reg '-MIPDECOconv-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
%            matlab2tikz(fullfile(figDir,[dataset '-' reg '-MIPDECOconv-' num2str(k) '-' num2str(j) '.tex']),'width','\iwidth','height','\iheight');   
%            return
           close(fig);
       end
   end
end

