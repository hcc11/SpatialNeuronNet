function MakeFigure5(varargin) 
% run Simulation_Fig5.m first 

addpath(genpath([pwd '/figtools/']));
%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',1); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = ''; 
setPlotOpt('plos','path',cDir,'cols',1,'height',15); 
inpath=[cDir 'data/']; 
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
% DC = axesDivide(2,[1 1 1 3],[0.11 0.08 0.82 0.85], 0.4, 0.5)';
DC = axesDivide(1,3,[0.17 0.07 0.7 0.90], .5, 0.5)';
% temp=axesDivide([1 1 1],2,[0.1 0.1 0.85 0.8], 0.4, 0.4)';
% DC(3:6)=temp(3:6);

Labels = {'A','B','C','D','E','F'}; LdPos = [-0.07,0.04];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
%       FigLabel(Labels{i},LdPos); 
end

fname1=[inpath 'neuralfield_stability'];
fname2=[inpath 'neuralfield_Lambda'];

HF_setFigProps;

% START PLOTTING 
data1=load(fname1);
data2=load(fname2);

[idy, idx]=find(data1.res(2).maxReWavNum<-0.1);
bdy_x=unique(idx);
bdy_y=zeros(2,numel(bdy_x));
for bb=1:numel(bdy_x)
    temp=find(idx==bdy_x(bb));
    bdy_y(1,bb)=min(idy(temp));
    bdy_y(2,bb)=max(idy(temp));
end

iA=1;
axes(AH(iA));
imagesc(data1.taui_use,data1.sigmai_use,data1.res(1).maxReWavNum)
xlabel('\tau_i')
ylabel('\sigma_i')
% axis square
wavenums=unique(data1.res(1).maxReWavNum); 
maxval = max(wavenums); %find maximum intensity
map = colormap; %get current colormap (usually this will be the default one)
imdata=data1.res(1).maxReWavNum;
ind=find(imdata<-0.1);
imdata = floor((imdata./maxval)*(length(map)-1))+1; 
imrgb=ind2rgb(imdata, map);
[x, y] = ind2sub(size(imrgb),ind);
% now set colors to grey
for indx = 1:length(x)
      imrgb(x(indx),y(indx),:) = [.5 .5 .5];
end
caxis([0 maxval])
h = narrow_colorbar('vert');
Pos = get(h,'Position');
set(h,'Position',Pos+[0.01,0,0,0],'YAxisLocation','right');
text(1.2,0.5,'Wavenumber w/ max Re(\lambda)','Units','n','Rotation',90,'Horiz','c',AxisLabelOpt{:});
set(h,'ytick',[0 1 2 3])

image(data1.taui_use,data1.sigmai_use,imrgb) 
axis xy
dx=data1.taui_use(2)-data1.taui_use(1); 
dy=data1.sigmai_use(2)-data1.sigmai_use(1);
xlim([data1.taui_use(1)-dx/2 data1.taui_use(end)+dx/2])
ylim([data1.sigmai_use(1)-dy/2 data1.sigmai_use(end)+dy/2])

plot(data1.taui_use(bdy_x),data1.sigmai_use(bdy_y(1,:)),'--k')
plot(data1.taui_use(bdy_x),data1.sigmai_use(bdy_y(2,:)),'--k')
plot(data1.taui_use(bdy_x(1))*[1 1],data1.sigmai_use(bdy_y(:,1)),'--k')
plot(data1.taui_use(bdy_x(end))*[1 1],data1.sigmai_use(bdy_y(:,end)),'--k')


Titles={'\sigma_i = \sigma_e','\sigma_i > \sigma_e'};
for ss=1:2
iA=ss+1;
axes(AH(iA));
colororder=copper(numel(data2.taui_use)); 
for j = 1:1:numel(data2.taui_use)
plot(data2.wavenums, data2.Lambdas(ss).maxreallam(j,:),'color',colororder(j,:))
text(.7,.1*(j-1)+.1,sprintf('%.1f ms',data2.taui_use(j)),'unit','n','color',colororder(j,:))
end
xlim([0 5])
plot(xlim,[0 0],'k--')
xlabel('Wavenumber')
ylabel('Max Re(\lambda)')
title(Titles{ss})
end

HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);

