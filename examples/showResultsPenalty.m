clear;

load('Peaks2D.mat')
dataset = 'Peaks2D';
reg     = 'wTVReg';
m       = [256 128];
domain = [0 2 0 1];
doPrint = 1;

files = dir([dataset '-' reg '-noise-0.1-' num2str(m(1)) 'x' num2str(m(2)) '-penalty.mat']);
figDir = '/Users/lruthot/Dropbox/Projects/2019-MIPDECOpaper/images/figPeaks2D/';

fname = files(1).name;
load(fname)


%% plot solution of SourcesMIPDECO
xa = linspace(domain(1),domain(2),size(W,1));
ya = linspace(domain(3),domain(4),size(W,2));
[XA,YA] = meshgrid(xa,ya);
xfine =getCellCenteredGrid(domain,size(W));
IoU = zeros(size(SourcesMIPDECOPenalty,3),size(SourcesMIPDECOPenalty,4),2);

for k=1:size(SourcesMIPDECOPenalty,3)    
    for j=1:size(SourcesMIPDECOPenalty,4)

        fig = figure(); clf;
        viewImage2Dsc(SourcesMIPDECOPenalty(:,:,k,j),domain,m);
        hold on;
        contour(XA',YA',W,1,'-r','LineWidth',3);
        caxis([0 1]);
        axis off
        
        mip = reshape(nnInter(SourcesMIPDECOPenalty(:,:,k,j),domain,xfine),size(W));
        intersect = mip.*W;
        union     = (mip+W)>eps;
        IoU(k,j,1) = nnz(intersect)/nnz(union);
        
        
        text(1.2,.1,sprintf('IoU=%1.1f%%',IoU(k,j,1)*100),'FontSize',40,'Color',[.99 .99 .99])
        if doPrint
            
            printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-SourcesMIPDECOPenalty-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
            close(fig);
        end
    end
end

%% plot initiliztion of SourcesMIPDECO
xa = linspace(domain(1),domain(2),size(W,1));
ya = linspace(domain(3),domain(4),size(W,2));
[XA,YA] = meshgrid(xa,ya);

for k=1:size(SourcesMIPDECOPenalty,3)
   for j=1:size(SourcesMIPDECOPenalty,4)

        fig = figure(); clf;
        viewImage2Dsc(InitMIPDECOPenalty(:,:,k,j),domain,m);
        hold on;
        contour(XA',YA',W,1,'-r','LineWidth',3);
        caxis([0 1]);
        
        mip = reshape(nnInter(InitMIPDECOPenalty(:,:,k,j),domain,xfine),size(W));
        intersect = mip.*W;
        union     = (mip+W)>eps;
        IoU(k,j,2) = nnz(intersect)/nnz(union);
        
        
        text(1.2,.1,sprintf('IoU=%1.1f%%',IoU(k,j,2)*100),'FontSize',40,'Color',[.99 .99 .99])
        
        
        axis off
        if doPrint
            printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-InitMIPDECOPenalty-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
            close(fig);
        end
   end
end
%%
figure;
plot(100*IoU(:,1,1),'linewidth',2,'DisplayName','TR');
hold on;
plot(100*IoU(:,2,1),'linewidth',2,'DisplayName','TR,local');
hold on;
plot(100*IoU(:,1,2),'linewidth',2,'DisplayName','penalty (init)');
lg = legend();
axis tight
lg.Location = 'SouthEast';
xlabel('iterations of penalty method')
ylabel('IoU')
set(gca,'FontSize',20)
if doPrint
	matlab2tikz(fullfile(figDir,[dataset '-' reg '-MIPDECOPenaltyIoU.tex']),'width','\iwidth','height','\iheight');   
end

%%%
%if doPrint
%   for k=1:size(HisMIPDECO,3)
%       for j=1:size(HisMIPDECO,4)
%           his = HisMIPDECO(:,:,k,j)
%           fig = figure(); clf;
%           iter = numel(find(his(:,1)));
%           yyaxis left
%           plot(his(1:iter,1),'linewidth',3)
%           ylabel('objective function')
%           
%           yyaxis right
%           plot(his(1:iter,3),'linewidth',3)
%           ylabel('Trust Region radius')
%           xlabel('iterations')
%           set(gca,'FontSize',20)
%           axis([0 50 0 512]);
%           printFigure(fig,fullfile(figDir,[dataset '-' reg '-MIPDECOconv-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
%%            matlab2tikz(fullfile(figDir,[dataset '-' reg '-MIPDECOconv-' num2str(k) '-' num2str(j) '.tex']),'width','\iwidth','height','\iheight');   
%%            return
%           close(fig);
%       end
%   end
%end

