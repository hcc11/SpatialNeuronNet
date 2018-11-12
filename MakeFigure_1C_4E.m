function MakeFigure_1C_4E(varargin) 
% Fig. 1C, Fig 4E & Fig. S6A-C,F-H
% no distance 
addpath(genpath([pwd '/figtools/']));

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',1); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = '';
setPlotOpt('plos','path',cDir,'cols',1.5,'height',16); 
inpath=[cDir 'data/']; 
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(3 ,4,[0.12 0.08 0.8 0.86], [0.4 .6], .4);
for kk=1:4
    iA=kk+4; 
    DC{iA}(4)=DC{iA}(4)*.75;
end
Labels = {'A1','B1','C1','D1','A2','B2','C2','D2','A3','B3','C3','D3'}; LdPos = [-0.07,0.04];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
%       FigLabel(Labels{i},LdPos); 
end
fname=[inpath 'GPFA_sim2'];
if P.Recompute  
    LF_generateData(fname,inpath);
    set(gca,'linewidth',1)
end 

HF_setFigProps;

% START PLOTTING 
% load(fname);
data(1)=load([inpath 'FA_data'],'Lambda','COVm','Qm','eigvector1');
data(2)=load([inpath 'FA_model_Jex25'],'Lambda','COVm','Qm','eigvector1');
data(3)=load([inpath 'FA_model_fastInh_Jex25'],'Lambda','COVm','Qm','eigvector1');
data(4)=load([inpath 'FA_model_broadInh_Jex25'],'Lambda','COVm','Qm','eigvector1');

eigvector{1}=[[data(1).eigvector1{:,1}]',[data(1).eigvector1{:,2}]'];
eigvector{2}=reshape(data(2).eigvector1(:,:,[3 6]),[],2);
eigvector{3}=reshape(data(3).eigvector1,[],2);
eigvector{4}=reshape(data(4).eigvector1,[],2);
data(2).Lambda=data(2).Lambda(:,:,[3 6]); 
data(2).COVm=data(2).COVm(:,[3 6]); 
data(2).Qm=data(2).Qm(:,[3 6]); 

colorAU= [0    0.5000    0.4000;
    0.9290    0.6940    0.1250];
Titles={'V4 data','model','model w/ fast Inh','model w/ broad Inh'};
YLims2={[0 20],[0 40],[0 10],[0 100]};
for kk=1:4
    iA=kk;
    axes(AH(iA));
    M=5;
    Ntrial=size(data(kk).Lambda,2);
    for pid=1:2
    errorbar(1:M,mean(data(kk).Lambda(1:M,:,pid),2),std(data(kk).Lambda(1:M,:,pid),[],2)/sqrt(Ntrial),'color',colorAU(pid,:),'linewidth',1)
    end
    ylim(YLims2{kk})
    xlim([0 5])
    set(gca,'xtick',0:5)
    set(gca,'ytick',[0 YLims2{kk}(2)/2 YLims2{kk}(2)])
    if kk==4
    xlabel('eigenmode')
    end
    ylabel('eigenvalue')
    title(Titles{kk})
    if kk==1
        text(.6,.9,'Unattended','unit','n','Horiz','left','color',colorAU(1,:),'fontsize',8)
        text(.6,.7,'Attended','unit','n','Horiz','left','color',colorAU(2,:),'fontsize',8)
    end
end

for kk=1:4
    iA=kk+4; 
    axes(AH(iA));
    edges=-.5:.05:.5;
    for pid=1:2
    [h, edges0]=histcounts(eigvector{kk}(:,pid),edges);
    stairs(edges(1:end-1),h/sum(h),'color',colorAU(pid,:),'linewidth',1)
    end
    if kk==1
    ylim([0 .3])
    else
        ylim([0 .24])
    end
    plot([0 0],ylim,'k--','linewidth',1)
    xlim([-.5 .5])
    set(gca,'xtick',[-.5 0 .5])
    set(gca,'ytick',[])
    title('1st mode')
    xlabel('projection weight') 
end   

YLims1={[-0.02 1],[0 .8],[0 .08],[-.2 .6]};
for kk=1:4
    iA=kk+8;
    axes(AH(iA));
    Ns=size(data(kk).COVm,1);
    dw=0.2; 
    Xoffset=[-.2 .2];
    for pid=1:2
        Xpos=[1 2]+Xoffset(pid);
        plot(ones(Ns,1)*Xpos+dw*(rand(Ns,2)-0.5), [data(kk).COVm(:,pid), data(kk).Qm(:,pid)],'.','color',colorAU(pid,:),'markersize',5)
        plot(Xpos(1)+[-dw  dw],mean(data(kk).COVm(:,pid))*[1 1],'k-','linewidth',1)
        plot(Xpos(2)+[-dw  dw],mean(data(kk).Qm(:,pid))*[1 1],'k-','linewidth',1)
          errorbar(Xpos(1),mean(data(kk).COVm(:,pid)),std(data(kk).COVm(:,pid)),'k-','CapSize',5)
        errorbar(Xpos(2),mean(data(kk).Qm(:,pid)),std(data(kk).Qm(:,pid)),'k-','CapSize',5)
    end
    ylim(YLims1{kk})
    xlim([0.5 2.5])
    ylabel('covariance')
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'raw','residual'})
end



HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);
