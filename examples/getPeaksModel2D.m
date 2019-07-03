clc; clear;

S = flipud(peaks(256))';
S = S>2;
S(20+(1:75),200+(1:75)) = S(75+(1:75),125+(1:75));
S(75+(1:75),125+(1:75)) = 0;
S = [S 0*S];
S = flipud(S)';


domain = [0 2 0 1];
m = size(S);
W = S; 
figure(1); clf;
imagesc(S);
axis equal

 save('Peaks2D.mat','W','domain','m')

%%
% load Peaks2D
x1 = linspace(domain(1),domain(2),m(1));
x2 = linspace(domain(3),domain(4),m(2));

figure(1); clf;
imagesc(x1,x2,flipud(W'))
axis equal tight
figDir = '/Users/lruthot/Dropbox/Projects/2018-MIPDECOpaper/images/2Dpeaks';
set(gca,'FontSize',20)
 printFigure(gcf,fullfile(figDir,'2Dpeaks-W.png'),'printFormat','.png','printOpts','-dpng');

 %%
x1 = linspace(domain(1),domain(2),m(1)+1);
x2 = linspace(domain(3),domain(4),m(2)+1);
 figure(2); clf;
imagesc(x1,x2,flipud(utrue'))
hold on;
plot(rec(:,1),rec(:,2),'.r','MarkerSize',20)
axis equal tight
figDir = '/Users/lruthot/Dropbox/Projects/2018-MIPDECOpaper/images/2Dpeaks';
set(gca,'FontSize',20)
 printFigure(gcf,fullfile(figDir,'2Dpeaks-utrue.png'),'printFormat','.png','printOpts','-dpng');
