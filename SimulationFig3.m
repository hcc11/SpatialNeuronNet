%  simulations for Fig. 3
%  two-layer network, for different spatial structure and tau_i
%  each simulation takes about 1.5 hours  

clear

Wtype='broadRec';  % Fig. 3Aii, Aiv
% Wtype='uniformW';  % Fig. 3Ai, Aiii

data_folder='data/';
dim='2D';
taui_range = [.5 2:1:15]; % time scale of inhibitory current (Fig. 2d)

%%%%%%%% run on cluster %%%%%%%%%%%
rng('shuffle');
AI = getenv('PBS_ARRAYID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
% job_dex range from 1 to Np 

Np=length(taui_range); 
taudsyni=taui_range(job_dex); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

opt.save=1; % save data 
opt.CompCorr=0; % compute correlations 
     Nc=[500 500]; 
     % # of E neurons sampled from Layer 2 & 1,  when opt.Layer1only=1
     % # of E neurons sampled from Layer 3 & 2,  when opt.Layer1only=0
opt.Layer1only=1; % 1 for two-layer network, 0 for three-layer network  
opt.loadS1=0;
opt.plotPopR=0; % plot population rate
opt.fixW=0; 
%     Wseed1=Wseed1_range(nws);
%     Wseed2=Wseed2_range(nws); 


dt=.01;
T=20000;  % total simulation time (ms) 
filename=strrep(sprintf('%sRF%s2layer_%s_tausyni%.03g',...
            data_folder,dim, Wtype, taudsyni),'.','d'),

% parameters to change, default parameters are in RF2D3layer.m
ParamChange={'filename', filename;...
    'dt', dt; 'T', T; 'Nc',Nc;'param(1).taudsyn(3)', taudsyni};
if strcmp(Wtype,'uniformW')
    ParamChange=[ParamChange;{'param(1).sigmaRX', 10*ones(2,1);
'param(1).sigmaRR',10*ones(2,2)}];
end
if opt.fixW
    ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
end

%%%%%%%% run simulation %%%%%%%%%%%
RF2D3layer(opt, ParamChange)

%% Compute average Correlation
load(filename)
tauf_range=[5 10 20 40]; % width of Gaussian filters (ms)
Nc=1000; % # of neurons to sample 
Ne1=sqrt(param(1).Ne);

I1=transpose(unique(s1(2,:)));
Ix=(ceil(I1/Ne1))/Ne1; % x,y indexes of neural location 
Iy=(mod((I1-1),Ne1)+1)/Ne1;
I1=I1(Ix<0.75 & Ix>0.25 & Iy<0.75 & Iy>0.25);
I1=randsample(I1,Nc);
% compute spike counts using sliding window
Tw=200;
Tburn=1000;
time=0:1:T; 
re0=zeros(Nc,length(time));
for mm=1:Nc
    re0(mm,:)=hist(s1(1,I1(mm)-1/4<s1(2,:) & s1(2,:)<=I1(mm)+1/4),time)*1e3;
end
t=0:400; filter='Gaussian';
for ff=1:length(tauf_range)
    tauf=tauf_range(ff);
    h=exp(-(t-200).^2/(2*tauf^2))/sqrt(2*pi)/tauf;  % filter  
    re1=imfilter(re0(:,Tburn+1:end),h);re1=re1(:,Tw/2-2:end-Tw/2); 
    ind1=mean(re1,2)>0;
    re1=re1(ind1,:);
    COV=cov(re1');
    R = corrcov(COV);
    U=triu(ones(size(R)),1);
    Cbar_tauf(ff)=mean(R(U==1)); % average corr.
    COVbar_tauf(ff)=mean(COV(U==1)); % average cov.
end
save(filename,'Cbar_tauf','COVbar_tauf','tauf_range','filter','Nc','-append')      

%%  Collect data after all simulations 

% taui_range = [.5 2:1:15]; 
% filename=@(Wtype,taudsyni) strrep(sprintf('data/RF%s2layer_%s_tausyni%.03g',...
%     dim, Wtype, taudsyni),'.','d');
% 
% Wtypes={'broadRec', 'uniformW'};
% 
% Ntau=length(taudi_range);
% for type=1:2
%     res(type).filename=cell(1,Ntau);
%     res(type).Corr=zeros(Ntau,1);
%     res(type).COV=zeros(Ntau,1);
%     res(type).nu=zeros(Ntau,2);
%     res(type).FF=zeros(Ntau,1);
%     res(type).re=cell(1,Ntau);
%     res(type).Cee_d=zeros(Ntau,20);
%     res(type).COVee_d=zeros(Ntau,20);
%     res(type).Cbar_tauf=zeros(Ntau,4);
%     res(type).COVbar_tauf=zeros(Ntau,4);
%     
%     for k=1:Ntau
%         taudsyni=taudi_range(k);
%         load(filename(Wtypes{type},taudsyni))
% 
%         time=0:1:T;Ne=200^2;
%         res(type).filename{k}=filename(Wtypes{type},taudsyni);
%         res(type).Corr(k)=Cbar(1);
%         res(type).COV(k)=COVbar(1);
%         res(type).nu(k,:)=nuSim(1:2);
%         res(type).FF(k)=mean(var2./rate2)/5;
%         res(type).re{k}=hist(s0(1,s0(2,:)<=Ne),time)/Ne*1e3;
%         res(type).Cee_d(k,:)=C(:,1)';
%         res(type).COVee_d(k,:)=COV_d(:,1)';
%         res(type).Cbar_tauf(k,:)=Cbar_tauf;
%         res(type).COVbar_tauf(k,:)=COVbar_tauf;
%     end
% end
% save('data/RF2D_broadRec_uniformW_tausyni_sum.mat','res','taudi_range','daxis','tauf_range','filter','Nc')
% 


