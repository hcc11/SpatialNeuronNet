% simulation of the three-layer network (Fig. 4) 
% vary the static input (muE, muI) to Layer 3 E & I neurons 
% meaure noise correlations within Layer 3, and those between Layer 2 and 3 
% each simulation takes about 3.5 hours and 3 gb memory  

clear

data_folder='data/';

testp.inE=[0.];
testp.inI=[0.1 .15 .2 .25 .3 .35 .4]; 

%%%%%%%%%%%% for job array on cluster %%%%%%%%%%%%%%%% 
rng('shuffle');
AI = getenv('PBS_ARRAYID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
% job_dex range from 1 to 350  % Ntrial=50 per inI 

Np=length(testp.inI)*length(testp.inE); 
pid=mod(job_dex-1,Np)+1;
trial=ceil(job_dex/Np); 
ip1=1;  
ip2=pid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

dim='2D';

opt.save=1; % save data 
opt.CompCorr=1; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 3 & 2  
opt.loadS1=0;
opt.plotPopR=0; % plot population rate
opt.fixW=0; 
%     Wseed1=Wseed1_range(nws);
%     Wseed2=Wseed2_range(nws); 

inE=testp.inE(ip1); % input to Layer 3 exc. neurons 
inI=testp.inI(ip2); % input to Layer 3 inh. neurons 
Jx=25*[1;0.6]; % ffwd strength from Layer 2 to 3
Prx=[ .05; .05];
dt=.01;   % (ms)
T=20000;  % (ms)
filename=sprintf('%sRF%s3layer_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g',...
    data_folder,dim,Jx(1),Jx(2),inE,inI,trial,dt);
filename=strrep(filename,'.','d')

ParamChange={'param(2).Iapp', [inE;inI]; 'filename', filename;...
    'dt', dt; 'T', T; 'Nc',Nc;'param(2).Jx',Jx; 'param(2).Prx',Prx};
if opt.loadS1
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end
if opt.fixW
    ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
end

RF2D3layer(opt, ParamChange)

%% run after finishing all simulations 
% data_folder='data/';
% dim='2D';
% dt=0.01;
% 
% Jx=25*[1;0.6];
% testp.inE=[0];
% testp.inI=[0.1 .15 .2 .25 .3 .35 .4];
% 
% Ntrial=50;
% Np=length(testp.inI)*length(testp.inE);
% 
% res.nu=zeros(Np,5,Ntrial);
% res.Corr=zeros(Np,3,Ntrial);
% res.COV=zeros(Np,3,Ntrial);
% res.FF=zeros(Np,1,Ntrial);
% res.Cee_d=zeros(Np,20,Ntrial);
% res.Cex_d=zeros(Np,20,Ntrial);
% res.COVee_d=zeros(Np,20,Ntrial);
% res.COVex_d=zeros(Np,20,Ntrial);
% res.inI=zeros(Np,1);
% res.inE=zeros(Np,1);
% res.Imean=zeros(Np,2,Ntrial);
% fnames=cell(Np,Ntrial); 
% 
% % T=2e4;
% % time=0:1:T;
% % Ne2=200^2;
% 
% fnamesave=sprintf('%sRF%s3layer_Jex%.03g_Jix%.03g_dt%.03g_inI_sum',...
%         data_folder,dim,Jx(1),Jx(2),dt); 
% fnamesave=strrep(fnamesave,'.','d');
% 
% for pid=1:Np
%     inE=testp.inE(1);
%     inI=testp.inI(pid);
%     res.inI(pid)=inI; 
%     res.inE(pid)=inE; 
%     for trial=1:Ntrial
%         filename=sprintf('%sRF%s3layer_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g',...
%             data_folder,dim,Jx(1),Jx(2),inE,inI,trial,dt);
%         filename=strrep(filename,'.','d'),
%         data=load(filename,'COVbar','Cbar','var2','rate2','nuSim','C','COV_d','daxis','param');
%         fnames{pid,trial}=filename; 
%         res.COV(pid,:,trial)=data.COVbar;
%         res.Corr(pid,:,trial)=data.Cbar;
%         res.FF(pid,:,trial)=mean(data.var2./data.rate2)/5;
%         res.nu(pid,:,trial)=data.nuSim;
%         res.Cee_d(pid,:,trial)=data.C(:,1)';
%         res.Cex_d(pid,:,trial)=data.C(:,2)';
%         res.COVee_d(pid,:,trial)=data.COV_d(:,1)';
%         res.COVex_d(pid,:,trial)=data.COV_d(:,2)'; 
%         res.Imean(pid,:, trial)=data.param(2).Imean; 
%     %     re{pid,trial}=hist(data.s2(1,data.s2(2,:)<=Ne2),time)/Ne2*1e3;
%     end
% end
% daxis=data.daxis;
% save(fnamesave,'res','testp','daxis','fnames','-append')

