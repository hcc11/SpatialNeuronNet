function MakeFigure7(varargin) 

addpath(genpath([pwd '/figtools/']));

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',1); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',1);  

% SETUP BASICS
cDir = '';
setPlotOpt('custom','path',cDir,'cols',2,'height',10); 
inpath=[cDir 'data/']; 
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide([2 1 2 1 1.5],4,[0.07 0.1 0.9 0.83], [.15 .8 .15 .6], 0.4)';
temp = axesDivide([2 1 2 1 1.5],[1 1 .8 1.2],[0.07 0.1 0.9 0.83], [.15 .8 .15 1], [0.4 .5 .6])';
DC=DC(:);
DC(5:5:20)=temp(5:5:20);
DC{15}=DC{10};
DC{15}(2)=DC{15}(2)+DC{15}(4)*.5;
DC{15}(1)=DC{15}(1)+DC{15}(3)*.5;
DC{15}(4)=DC{15}(4)*.5;
DC{15}(3)=DC{15}(3)*.5;

DC{20}(4)=DC{20}(4)*.9;
DC{21}=DC{20};
DC{20}(2)=DC{20}(2)+DC{20}(4)*1.1;

Labels = {'A','B','C','D','E','F'}; LdPos = [-0.065,0.02];

for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
%       FigLabel(Labels{i},LdPos); 
%       AH2(i) = axes('Pos',DC2{i}); hold on; 
end

T0=1864+1.1e4;

fname1=[inpath 'RF2D3layer_fixW_Jex25_Jix15_Ntrial10_Nrep20_Re.mat'];
fname2=[inpath 'RF2D3layer_fixW_Jex25_Jix15_MacroChaos_rasters.mat'];

if P.Recompute  
    LF_generateData(fname1);
    LF_ffwdinput(fname2,T0)
end 

HF_setFigProps;

Titles={'Unattended','Attended'};
colorAU= [0    0.5000    0.4000;
    0.9290    0.6940    0.1250];

% START PLOTTING
load(fname1)
load(fname2)

Ne1=200;
Nrep=20;
tf=-100:100;
sig_f=2;
hf=exp(-tf.^2/(2*sig_f^2));
hf=hf/sum(hf);

for pid=1:2
    for row=1:4
        if row<=3
            iA=(row-1)*5+pid*2-1;
            axes(AH(iA))
            re2=imfilter(Re2{pid,Reps(row),trial},hf);
            re2=re2(t1:t2);
            plot((t1:t2)-t1,re2)
            xlim([0 t2-t1])
            ylim([0 150])
            set(gca,'xtick',[])
            plot((T0-t1)*[1 1],ylim,'k--')
            if row==1
                title(Titles{pid})
            end
            if pid==1 && row==2    
                h_text=text(-.3,.5,'popuplation rate (sp/s)','unit','n','color','k','horiz','center');
                set(h_text, 'rotation', 90)
            end
            if pid==1
                text(0.2,0.8,sprintf('trial %d',row),'unit','n','color','k','horiz','center');
            end
            
            iA=(row-1)*5+pid*2;
            axes(AH(iA))
            s2=res(pid,row).s2;
            Is=find(s2(1,:)<=T0 & s2(1,:)>T0-2 & s2(2,:)<Ne1^2);
            x=ceil(s2(2,Is)/Ne1);
            y=mod(s2(2,Is)-1,Ne1)+1;
            plot(x,y,'.','color',colorAU(pid,:),'markersize',2)
            axis([0 200 0 200])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            ht=title(sprintf('%.0f ms',T0-t1));
            axis square
            box on
        else
            iA=(row-1)*5+pid*2-1;
            axes(AH(iA))
            p = randperm(Nrep);
            for nrep=1:Nrep
                re2=imfilter(Re2{pid,p(nrep),trial}-Re2m{pid,trial},hf);
                re2=re2(t1:t2);
                plot(((t1:t2)-t1)*1e-3,re2)
                xlim([0 t2-t1]*1e-3)
            end
            ylim([-100 200])
            set(gca,'ytick',[-100 0 100 200])
            if pid==1
                h_text=text(-.3,.5,'difference (sp/s)','unit','n','color','k','horiz','center');
                set(h_text, 'rotation', 90)
            end
            xlabel('time (sec)')
            
            iA=(row-1)*5+pid*2;
            axes(AH(iA))
            Jx=25/sqrt(5e4); 
            imagesc(Input*Jx)
            axis xy
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            axis([.5 200.5 0 200.5])
            colormap('hot') 
            caxis([3 6])
            h = narrow_colorbar('vert');
            Pos = get(h,'Position');
            set(h,'Position',Pos+[0.0,0,0,0],'YAxisLocation','right');
            if pid==1
            text(1.4,0.5,'(pA)','Units','n','Rotation',90,'Horiz','c',AxisLabelOpt{:})
            end
            title('ffwd input')
            axis square
            box on
        end
    end
end

T=20; 
Ntrial=size(Rmean_s,2); 
for pid=1:2
    iA=pid*5;
    axes(AH(iA));
    for trial=1:Ntrial 
    Rvar_filter=imfilter(Rvar{pid,trial},hf');
    plot((1:length(Rvar{pid,trial}))*Nstep*1e-3+1+T*(trial-1),Rvar_filter,'color',colorAU(pid,:))
    end
    xlim([0 T*Ntrial])
    ylim([0 50])
    set(gca,'ytick',[0 50])
    if pid==1
        set(gca,'xtick',[])
        title('Var(E|X)')
    else
        plot([10,30],[20 20],'k')
        plot([50,50],[20 30],'k')
    end
    axis off
end

iA=3*5;
axes(AH(iA));
for pid=1:2
    mvar=0;
    for trial=1:Ntrial
        mvar=mvar+mean(Rvar{pid,trial});
    end
    mvar=mvar/Ntrial;
    hb=bar(pid,mvar,'Facecolor',colorAU(pid,:));
end
set(gca,'xtick',[])
xlim([0 3])
ylim([0 4])
set(gca,'ytick',[0 4])

iA=4*5;
axes(AH(iA));
XLims=[15 24];
YLims=[32 45];
Ntrial=size(Rmean_s,2); 
for pid=[1 2]
    for trial=1:Ntrial
    axes(AH(iA+pid-1));
    plot(re1syn_s{pid,trial},Rmean_s{pid,trial},'.','color',colorAU(pid,:),'markersize',2)
    end
axis([XLims YLims])
set(gca,'ytick',YLims)
if pid==1
    set(gca,'xtick',[]);
else
    set(gca,'xtick',XLims)
end
end
xlabel('X (sp/s)')
text(-.3,1.2,'R (sp/s)','unit','n','color','k','horiz','center', 'rotation', 90);

HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);


function  LF_ffwdinput(fname,T0)
load(fname)
s1=s1(:,s1(1,:)<T0&s1(1,:)>(T0-500));
Input_s1=zeros(200,200);
t=linspace(0,500,5001);
epsc= @(t) 0.2*(exp(-t./param(2).taudsyn(1,1))-exp(-t./param(2).taursyn(1,1)))/(param(2).taudsyn(1,1)-param(2).taursyn(1,1))+0.8*(exp(-t./param(2).taudsyn(1,2))-exp(-t./param(2).taursyn(1,2)))/(param(2).taudsyn(1,2)-param(2).taursyn(1,2));
for k=1:200^2
tsps=s1(1,s1(2,:)<k+0.1&s1(2,:)>k-0.1);
Input_s1(k)=sum(epsc(T0-tsps));
end

pex0=param(2).Prx(1);
pix0=param(2).Prx(2);
Ne=200^2;Ni=100^2;
Kex=ceil(pex0*Ne);
Kix=ceil(pix0*Ni);
Kx=Kex+Kix;
rng(Wseed)
[Wrr2,Wrf2]=gen_weights(param(2).Ne,param(2).Ni,param(2).Nx,param(2).sigmaRX,param(2).sigmaRR,param(2).Prr,param(2).Prx,'2D');
clear Wrr2
Input=zeros(200,200);
for k=1:200^2
post=Wrf2((k-1)*Kx+1:(k-1)*Kx+Kex);
[c,post]=hist(post,unique(post));
Input(post)=Input(post)+c'.*Input_s1(k);
end
rng('shuffle')
save(fname,'Input_s1','Input','-append')

function  LF_generateData(fname)
load(fname)
t_h=0:500;
hs=(exp(-t_h/100)-exp(-t_h/2))./(100-2);
hf=(exp(-t_h/5)-exp(-t_h/1))./(5-1);
Tw=200; % sliding window size
Tburn=1000;
Nstep=10; %step size for sliding window
Nt=ceil((T-Tburn-Tw)/Nstep);
Re2_s=zeros(Nt,Nrep); 
re1syn_s=cell(Np,Ntrial);
Rvar=cell(Np,Ntrial);
Rmean_s=cell(Np,Ntrial);
Re2m=cell(Np,Ntrial);
for pid=1:Np
    for trial=1:Ntrial
        Re2m{pid,trial}=zeros(size(Re2{pid,1,trial})); 
        for nrep=1:Nrep
            re2_s=(imfilter(Re2{pid,nrep,trial}(Tburn+1:end),ones(1,Tw)/Tw));
            re2_s=re2_s(Tw/2+1:Nstep:end-Tw/2-1);
            Re2_s(:,nrep)=re2_s';
            Re2m{pid,trial}=Re2m{pid,trial}+Re2{pid,nrep,trial};
        end
        Re2m{pid,trial}=Re2m{pid,trial}/Nrep;
        Rvar{pid,trial}=var(Re2_s,[],2);
        Rmean_s{pid,trial}=mean(Re2_s,2);
        re1syn_s{pid,trial}=Re1{pid,trial};
        re1syn_s{pid,trial}=imfilter(re1syn_s{pid,trial}(Tburn+1:end),ones(1,Tw)/Tw);
        re1syn_s{pid,trial}=re1syn_s{pid,trial}(Tw/2+1:Nstep:end-Tw/2-1);
    end
end
save(fname,'Rvar','Rmean_s','re1syn_s','Nstep','Re2m','-append')


