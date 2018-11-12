function MakeFigure3(varargin) 

addpath(genpath([pwd '/figtools/']));

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',1); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',1);  

% SETUP BASICS
cDir = ''; 
setPlotOpt('plos','path',cDir,'cols',2,'height',13); 
inpath=[cDir 'data/']; 
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide([0.6 1.5 1.5 1.],[.6 1.5 3],[0.1 0.08 0.85 0.9], [.3 .6 .6],[0.5 .7 ])';
temp = axesDivide([0.6 1.5 1.5 1.],[1 1 1 1 2],[0.1 0.08 0.85 0.85], [.3 .6 .6],[0.2 1 .2 1])';
DC3=temp(4,:);

DC([1 4],:)=[];

Labels = {'A','','D','B','','E','C','','F'}; LdPos = [-0.07,0.02];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
%       FigLabel(Labels{i},LdPos); 
end
for i = 1:numel(DC3)
    AH3(i) = axes('Pos',DC3{i}); hold on;
end
for iA=[5 6]
    axes(AH(iA));
    axis off
    DC2=DC(iA);
    DC2{1}(3)=DC2{1}(3)*.6;
    DC2{1}(1)=DC2{1}(1);
    DC2{1}(4)=DC2{1}(4)/2;
    DC2{1}(2)=DC2{1}(2)+DC2{1}(4);
    DC2{2}=DC2{1};
    DC2{2}(2)=DC2{1}(2)-DC2{1}(4)*.6;
    DC2{2}(1)=DC2{1}(1)+DC2{1}(3)*.4;
    DC2{3}=DC2{2};
    DC2{3}(2)=DC2{2}(2)-DC2{2}(4)*.6;
    DC2{3}(1)=DC2{2}(1)+DC2{2}(3)*.4;
    for i = 1:numel(DC2)
        AH2(i,iA) = axes('Pos',DC2{i}); hold on;
    end
end

HF_setFigProps;

% sample raster
Nc=200; % # of sampled neurons 
t1=1000;t2=1400;
if P.Recompute  
    LF_generateData([inpath 'RF2D_uniformW_tausyni8'],Nc,t1,t2);
    LF_generateData([inpath 'RF2D_uniformW_tausyni1'],Nc,t1,t2);
end

% START PLOTTING 

T0=[4600 4618 4628]; % snapshot time 

data{2,2}=load([inpath 'RF2D_broadRec_tausyni8'],'s0','rate','R');
data{2,1}=load([inpath 'RF2D_broadRec_tausyni1'],'s0','rate','R');
data{1,2}=load([inpath 'RF2D_uniformW_tausyni8'],'ts','Ic','rate','R','taursyni','taudsyni','taursyne','taudsyne');
data{1,1}=load([inpath 'RF2D_uniformW_tausyni1'],'ts','Ic','rate','R','taursyni','taudsyni','taursyne','taudsyne');

Ne1=200;
colororder=lines(3);

for jj=1:2
    iA=jj; 
    axes(AH(iA));
    taursyni=data{1,jj}.taursyni;taudsyni=data{1,jj}.taudsyni;
    taursyne=data{1,jj}.taursyni;taudsyne=data{1,jj}.taudsyne;
    syn_t=linspace(0,30,201);   
    synE=(exp(-syn_t./taursyne)-exp(-syn_t./taudsyne))/(taursyni-taudsyne);
    synI=(exp(-syn_t./taursyni)-exp(-syn_t./taudsyni))/(taursyni-taudsyni);
    plot(syn_t,synE,'-','color',colororder(jj,:));
    plot(syn_t,synI,'--','color',colororder(jj,:));
    if jj==1
        plot([15 20], [.5 .5],'k')
        text(23,.5,'EPSC','unit','data','color','k')
        plot([15 20], [.3 .3],'--k')
        text(23,.3,'IPSC','unit','data','color','k')
        plot([0 5],[-.1 -.1],'k')
        text(0,-0.3,'5 ms','unit','data','color','k')   
    end
    ylim([-.1 .5])
    axis off
end

ii=1;
for jj=1:2
    iA=2*ii+jj;
    axes(AH(iA));
    ts=data{ii,jj}.ts;
    Ic=data{ii,jj}.Ic;
    for mm=1:Nc 
        plot(ts{mm},mm*ones(size(ts{mm})),'k.','markersize',3)
    end
    axis([t1 t2 1 Nc])
    xlabel('time (ms)')
    ylabel('neuron ID')
end
ii=2;
for jj=1:2
    iA=2*ii+jj;
    s0=data{ii,jj}.s0;
    for k=1:3
        t1=T0(k);
        axes(AH2(k,iA))
        Is=find(s0(1,:)<=t1 & s0(1,:)>t1-1 & s0(2,:)<Ne1^2);
        x=ceil(s0(2,Is)/Ne1);
        y=mod(s0(2,Is)-1,Ne1)+1;
        plot(x,y,'k.','markersize',3)
        axis([0 200 0 200])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        if k==1
            ht=title(sprintf('%.0f ms',t1));
        else
            ht=title(sprintf('%.0f ms',t1));
            pos=get(ht,'Position');
            pos(1)=pos(1)+85;
            set(ht,'Position',pos)
        end
        box on
        axis square
        if k==1
            ylabel('neuron location (Y)')
        elseif k==3
            xlabel('neuron location (X)')
        end
    end
end

k=1; 
for ii=1:2
axes(AH3(2*(k-1)+ii));
[h,c]=hist([data{ii,1}.rate,data{ii,2}.rate],0:1:70);
h=h./(ones(size(h,1),1)*sum(h,1));
for jj=1:2
    plot(c,h(:,jj),'color',colororder(jj,:))
end
axis([0 70 0 .06])
set(gca,'ytick',[])
set(gca,'xtick',[0 70])
if ii==1
    set(gca,'xtick',[])
else
    xlabel('firing rate (sp/s)')
    text(-.18,1.1,'P(r)','unit','n','color','k','horiz','center','rotation', 90);
end
end

k=2;
for ii=1:2
    axes(AH3(2*(k-1)+ii));
    [h,c]=hist([data{ii,1}.R,data{ii,2}.R],-1:.02:1);
    h=h./(ones(size(h,1),1)*sum(h,1))/.02;
    for jj=1:2
        plot(c,h(:,jj),'color',colororder(jj,:))
        plot(mean(data{ii,jj}.R)*[1 1],[0 4.5],'--','color',colororder(jj,:))
    end
    axis([-1 1 0 4.5])
    set(gca,'ytick',[])
    if ii==1
        set(gca,'xtick',[])
    else
        xlabel('correlation')
        text(-.18,1.1,'P(corr)','unit','n','color','k','horiz','center','rotation', 90);
    end
end

axes(AH3(5));
data_taui=load([inpath 'RF2D_broadRec_uniformW_tausyni_sum.mat']);
plot(data_taui.taudi_range,data_taui.res(1).Cbar_tauf(:,2),'k')
plot(data_taui.taudi_range,data_taui.res(2).Cbar_tauf(:,2),'--k')
plot([7 10], [.8 .8],'k')
text(11,.8,'2D','unit','data','color','k')
plot([7 10], [.65 .65],'--k')
text(11,.65,'0D','unit','data','color','k')
xlabel('\tau_i')
text(-.18,.5,'correlation','unit','n','color','k','horiz','center','rotation', 90);
xlim([0 15])
set(gca,'ytick',[0 1])


HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);

function  LF_generateData(fname,Nc,t1,t2)
load(fname,'s0','Ne1')
Ic=randsample(Ne1^2,Nc);
ts=cell(size(Ic));
for mm=1:Nc
    ts{mm}=s0(1,(s0(1,:)<=t2 & s0(1,:)>t1 & Ic(mm)-1/4<s0(2,:) & s0(2,:)<=Ic(mm)+1/4));
end
save(fname,'ts','Ic','-append')

