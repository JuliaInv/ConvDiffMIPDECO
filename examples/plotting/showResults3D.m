clc;clear; 

load('Sources3D.mat')
dataset = 'Sources3D';
 reg     = 'wTVReg';
%  reg     = 'wTVReg';
m = [96 48 48];
 doPrint = 1;



files = dir([dataset '-' reg '-noise-0.1-*.mat']);
figDir = '/Users/lruthot/Dropbox/Projects/2019-MIPDECOpaper/images/figSource3D/';


 
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
src = reshape(Sources(:,idx),m);
viewSlices(src,domain,m,'s2',35,'s3',12,'s1',50);
caxis([0 1])
axis equal tight

%% 
if doPrint
    for k=1:size(Sources,2);
       src = reshape(Sources(:,k),m);

       fig = figure(); clf;
       viewSlices(src,domain,m,'s2',35,'s3',12,'s1',50);
       caxis([0 1])
       axis off
axis equal tight
       printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(k) '.png']),'printOpts','-dpng','printFormat','.png');
       close(fig);
    end
end

%% plot solution of SourcesMIPDECO
xa = linspace(domain(1),domain(2),size(W,1));
ya = linspace(domain(3),domain(4),size(W,2));
xfine = getCellCenteredGrid(domain,size(W));
[XA,YA] = meshgrid(xa,ya);
if doPrint
    for k=1:size(SourcesMIPDECO,4)
        for j=1:size(SourcesMIPDECO,5)
            fig = figure(); clf;
            mip = SourcesMIPDECO(:,:,:,k,j);
            viewSlices(mip,domain,m,'s2',35,'s3',12,'s1',50);
            mip = reshape(nnInter(mip,domain,xfine),size(W));
            intersect = mip.*W;
            union     = (mip+W)>eps;
            IoU = nnz(intersect)/nnz(union);
            
            axis equal off
            text(1.2,.2,sprintf('IoU=%1.1f%%',IoU*100),'FontSize',35,'Color','k')
            
            caxis([0 1]);
            axis off
            printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-SourcesMIPDECO-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
            close(fig);
        end
    end
end

%% plot initiliztion of SourcesMIPDECO
xa = linspace(domain(1),domain(2),size(W,1));
ya = linspace(domain(3),domain(4),size(W,2));
[XA,YA] = meshgrid(xa,ya);
if doPrint
    for k=1:size(SourcesMIPDECO,4)
        for j=1:size(SourcesMIPDECO,5)
            fig = figure(); clf;
            mip = InitMIPDECO(:,:,:,k,j);
            viewSlices(mip,domain,m,'s2',35,'s3',12,'s1',50);
            mip = reshape(nnInter(mip,domain,xfine),size(W));
            intersect = mip.*W;
            union     = (mip+W)>eps;
            IoU = nnz(intersect)/nnz(union);
            
            axis equal off
            text(1.2,.2,sprintf('IoU=%1.1f%%',IoU*100),'FontSize',35,'Color','k')
            
            printFigure(fig,fullfile(figDir,[dataset '-' reg '-' num2str(m(1)) 'x' num2str(m(2)) '-InitMIPDECO-' num2str(k) '-' num2str(j) '.png']),'printOpts','-dpng','printFormat','.png');
            close(fig);
        end
    end
end

%%
if doPrint
   for k=1:size(HisMIPDECO,3)
       for j=1:size(HisMIPDECO,4)
           His = HisMIPDECO(:,:,k,j);
           fig = figure(); clf;
           iter = numel(find(HisMIPDECO(:,1)));
           yyaxis left
           plot(HisMIPDECO(1:iter,1),'linewidth',3)
           ylabel('objective function')
           
           yyaxis right
           plot(HisMIPDECO(1:iter,3),'linewidth',3)
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

