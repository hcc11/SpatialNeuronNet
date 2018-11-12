function MakeFigure6(varargin) 

addpath(genpath('~/Dropbox/Matlab_codes/BE_figTools/'))
%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',1); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = '';
setPlotOpt('custom','path',cDir,'cols',1,'height',9); 
inpath=[cDir 'data/']; 
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(2,2,[0.13 0.13 0.8 0.8], 0.5, .6)';

Labels = {'A','B','C','D'}; LdPos = [-0.07,0.04];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
      set(gca,'linewidth',1)
      set(gca,'fontsize',8)
%       FigLabel(Labels{i},LdPos); 
end

if P.Recompute  
    LF_generateData(fname,inpath);
end

HF_setFigProps;

% START PLOTTING 
% load(fname);

data=load([inpath 'FA_model_slowInh_Jex25_dist']);
load([inpath 'FA_data_dist'])

colorAU= [0    0.5000    0.4000;
    0.9290    0.6940    0.1250];

M=5;
colororder=copper(M);

%%%%%%%%%%%%  V4 data  %%%%%%%%%%%%%%%%%%%
% discard data points at distances at zero and >2mm
d1=d1(2:14);  % distance between electrodes
cov_d=cov_d(2:14,:); 
cov_d_std=cov_d_std(2:14,:); 
cov_d_Npair=cov_d_Npair(2:14,:); 
LL_d=LL_d(2:14,:,:); 
LL_d_std=LL_d_std(2:14,:,:); 
LL_d_Npair=LL_d_Npair(2:14,:,:); 

% average across sessions 
for pid=1:2
    for kk=1:length(d1)
        for mm=1:M
            for trial=1:72
                LL_d_tot{kk,mm,trial,pid}=LL_d_tot{kk,mm,trial,pid}*Lambda(mm,trial,pid);
            end
            LL_d(kk,mm,pid)=mean([LL_d_tot{kk,mm,:,pid}]);
            LL_d_std(kk,mm,pid)=std([LL_d_tot{kk,mm,:,pid}]);
            LL_d_Npair(kk,mm,pid)=length([LL_d_tot{kk,mm,:,pid}]);
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iA=1;
axes(AH(iA));
M=5;
Ntrial=size(data.cov_d,2);
for pid=1:2
    shadedErrorBar(data.daxis,mean(data.cov_d(1:20,:,pid),2),std(data.cov_d(1:20,:,pid),[],2)./sqrt(Ntrial),{'color',colorAU(pid,:),'linewidth',1})
end
ylim([-.1 .8])
xlim([0 .5])
text(.6,.9,'Unattended','unit','n','Horiz','left','color',colorAU(1,:),'fontsize',8)
text(.6,.7,'Attended','unit','n','Horiz','left','color',colorAU(2,:),'fontsize',8)
set(gca,'xtick',[0,.25 .5])
xlabel('distance (a.u.)')
ylabel('Covariance')
title('Model')

iA=3;
axes(AH(iA));
for mm=1:M
    data.LL_d{1}(1:20,mm,:)=data.LL_d{1}(1:20,mm,:).*permute(repmat(data.Lambda(mm,:,1),[20,1]),[1,3,2]);
    shadedErrorBar(data.daxis,mean(data.LL_d{1}(1:20,mm,:),3),...
        std(data.LL_d{1}(1:20,mm,:),[],3)./sqrt(Ntrial),{'color',colororder(mm,:),'linewidth',1})  
    text(1.,1-.1*mm,sprintf('%d',mm),'unit','n','color',colororder(mm,:),'fontsize',8)
end
text(1.,1,'mode','Horiz','c','unit','n','color','k','fontsize',8)
xlim([0 .5])
ylim([-.05 0.2 ])
set(gca,'xtick',[0,.25 .5])
xlabel('distance (a.u.)')
ylabel('\lambda_iv_i*v_i^T')

iA=2;
axes(AH(iA));
for pid=1:2
shadedErrorBar(d1*.4,cov_d(:,pid),cov_d_std(:,pid)./sqrt(cov_d_Npair(:,pid)),{'color',colorAU(pid,:),'linewidth',1})
end
title('V4 data')
xlabel('distance (mm)')
set(gca','xtick',0:1:3)
xlim([0 2.1])
ylim([-0.05 .4]) 

iA=4;
axes(AH(iA));
for mm=1:M
shadedErrorBar(d1*.4,LL_d(:,mm,1),LL_d_std(:,mm,1)./sqrt(LL_d_Npair(:,mm,1)),{'color',colororder(mm,:),'linewidth',1})
end
xlim([0 2.1])
ylim([-.05 0.3])
set(gca','xtick',0:1:3)
xlabel('distance (mm)')
 

 






HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);
