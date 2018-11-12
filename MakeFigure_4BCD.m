function MakeFigure_4BCD(varargin) 

addpath(genpath([pwd '/figtools/']));

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',1); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = '';
setPlotOpt('manuscript','path',cDir,'cols',1.5,'height',12); 
inpath=[cDir 'data/']; 
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
pos=[0.1 0.1 0.85 0.85];
DC = axesDivide(2,2,pos, 0.4, .6)';
DC(1)=DC(2); 
DC{1}(4)=DC{1}(4)*.45;
DC{2}(4)=DC{2}(4)*.45;
DC{1}(2)=DC{1}(2)+DC{1}(4)/.45*.55;


Labels = {'A','B','C','D','E','F','G'}; LdPos = [-0.07,0.02];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
%       FigLabel(Labels{i},LdPos); 
end
HF_setFigProps;

colorAU= [0    0.5000    0.4000;
    0.9290    0.6940    0.1250];
state={'Unattended','Attended'};

fnames=[inpath 'Fig3data.mat'];
load(fnames) 
Ntrial=size(Corr,3);
Np=length(testp.inI); 

% START PLOTTING 

for iA=1:2 
axes(AH(iA));
plot(0:1e-3:20,Re{iA},'color',colorAU(iA,:))
axis([0 20 0 400])
text(.5,.9,state{iA},'unit','n','Horiz','center','color','k')
xlabel('Time (sec)')
if iA==2
text(-0.2,1.2,'Pop. rate (Hz)','unit','n','Horiz','center','color','k','Rot',90)
end
end

iA=3; 
axes(AH(iA));
errorbar(testp.inI,mean(Corr(:,1,:),3),std(Corr(:,1,:),[],3)/sqrt(Ntrial),'k','linewidth',1)
plot(testp.inI(3), mean(Corr(3,1,:),3),'o','color',colorAU(1,:))
plot(testp.inI(6), mean(Corr(6,1,:),3),'o','color',colorAU(2,:))
axis([.05 .45 .055 .1])
set(gca,'ytick',[.06 .08 .1])
xlabel('Attention drive to inh') 
ylabel('Correlation (MT-MT)') 

iA=4; 
axes(AH(iA));
errorbar(testp.inI,mean(Corr(:,2,:),3),std(Corr(:,2,:),[],3)/sqrt(Ntrial),'k','linewidth',1)
plot(testp.inI(3), mean(Corr(3,2,:),3),'o','color',colorAU(1,:))
plot(testp.inI(6), mean(Corr(6,2,:),3),'o','color',colorAU(2,:))
axis([.05 .45 .015 .04])
set(gca,'ytick',[.02 .03 .04])
xlabel('Attention drive to inh') 
ylabel('Correlation (V1-MT)')


HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);


