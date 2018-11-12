% for simulations to study macroscopic chaos (Fig. 7)
% fix the spike trains from Layer 2 (s1) and connectivities (Wrr2, Wrf2)

clear

data_folder='data/';

Wseed1_range=[8541; 45134; 48395; 3547; 14845;  71109; 99911; 98570;...
       68790; 16203 ];
   
Wseed2_range=[800281; 141887; 421762; 915736; 792208; 959493; ...
    157614;  970593; 957167; 485376]; 

testp.inE=[0];
testp.inI=[.2 .35]; 

%%%%%%%%%%%% for job array on cluster %%%%%%%%%%%%%%%% 
rng('shuffle');
AI = getenv('PBS_ARRAYID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
% job_dex from 1 to 400 (Nrep*Ntrial*Np*Nnws)

Np=length(testp.inI)*length(testp.inE); 
Nrep=20;  % # of repeats per condition 
Ntrial=10; % # of realizations of connectivity and Layer 2 spike trains 

pid=mod(job_dex-1,Np)+1;
trial=mod(ceil(job_dex/Np)-1,Ntrial)+1;
nrep=ceil(job_dex/(Np*Ntrial));
nws=trial; 
ip1=1;
ip2=pid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

dim='2D';

opt.save=1; % save data 
opt.CompCorr=0; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.loadS1=1; % use stored spike trains from Layer 2
opt.plotPopR=0; % plot population rate
opt.fixW=1; 
    Wseed1=Wseed1_range(nws);
    Wseed2=Wseed2_range(nws); 
opt.savecurrent=0; 

inE=testp.inE(ip1);
inI=testp.inI(ip2);
Jx=25*[1; 0.6];
Prx=[ .05; .05];
dt=.01;
T=20000;  % (ms)
s1_fname=strrep(sprintf('%sRF%s3layer_fixW_slowInh_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g_nws1',...
    data_folder,dim,Jx(1),Jx(2),inE,.2,trial,dt),'.','d');
    
filename=strrep(sprintf('%sRF%s3layer_fixW_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_dt%.03g_nws%.03g_nrep%.03g',...
    data_folder,dim,Jx(1),Jx(2),inE,inI,dt,nws,nrep),'.','d');

ParamChange={'param(2).Iapp', [inE;inI]; 'filename', filename;...
    'dt', dt; 'T', T; 'Nc',Nc;'param(2).Jx',Jx; 'param(2).Prx',Prx};
if opt.loadS1
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end
if opt.fixW
    ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
end

RF2D3layer(opt, ParamChange)

%% collect population rates after all simulations  
% data_folder='data/';
% dim='2D';
% inE=0;
% inI=.2;
% Jx=25*[1; 0.6];
% Prx=[ .05; .05];
% dt=0.01;
% Nrep=20;
% Ntrial=10;
% T=2e4;
% time=0:1:T;
% Ne2=200^2;
% Ne1=200^2;
% Np=1; 
% pid=1; 
% fnamesave=sprintf('%sRF2D3layer_fixW_Jex%.03g_Jix%.03g_Ntrial%.03g_Nrep%.03g_Re',...
%             data_folder,Jx(1),Jx(2),Ntrial,Nrep);
% fnamesave=strrep(fnamesave,'.','d')
% for trial=1:Ntrial
%     nws=trial;
%     for nrep=1:Nrep
%         filename=sprintf('%sRF%s3layer_fixW_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_dt%.03g_nws%.03g_nrep%.03g',...
%             data_folder,dim,Jx(1),Jx(2),inE,inI,dt,nws,nrep);
%         filename=strrep(filename,'.','d')
%         data=load(filename,'s2','p_stim'); 
%         Re2{pid,nrep,trial}=hist(data.s2(1,data.s2(2,:)<=Ne2),time)/Ne2*1e3;
%     end
%     s1_fname=data.p_stim.s1_fname;
%     load(s1_fname,'s1')
%     Re1{pid,trial}=hist(s1(1,s1(2,:)<=Ne1),time)/Ne1*1e3;
% end
% save(fnamesave,'Re2','Re1','Ntrial','Nrep','Np','Jx','Prx','Ne1','Ne2','inE','inI','T')

%% collect spike times (s2, s1) for raster plots (Fig. 7B,C)
% data_folder='data/';
% trial=7; % sim # of Layer 2
% Reps=[9 10 12]; % repetition # of Layer 3
% t1=1.1e4; % start time; 
% t2=1.4e4; % end time 
% 
% testp.inI=[.2 .35]; % Unatt: inI=.2; Att: inI=.35;
% Jx=25*[1; 0.6];
% inE=0;
% dt=0.01; 
% nws=trial; 
% Ne=200^2; 
% fnamesave=sprintf('%sRF2D3layer_fixW_Jex%.03g_Jix%.03g_MacroChaos_rasters',data_folder,Jx(1),Jx(2));
% for pid=1:2
%     inI=testp.inI(pid);
%     for kk=1:3
%         nrep=Reps(kk);
%         filename=sprintf('%sRF2D3layer_fastslowJx_fixW_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_dt%.03g_nws%.03g_nrep%.03g',...
%             data_folder,Jx(1),Jx(2),inE,inI,dt,nws,nrep);
%         filename=strrep(filename,'.','d')
%         load(filename,'s2','p_stim','param');
%         res(pid,kk).s2=s2(:,s2(1,:)>t1&s2(1,:)<t2&s2(2,:)<=Ne);
%         res(pid,kk).filename=filename;
%     end
% end
% Wseed=param(2).Wseed; 
% load(p_stim.s1_fname,'s1')
% s1=s1(:,s1(1,:)>t1&s1(1,:)<t2&s1(2,:)<=Ne); % s1 is the same for Unatt and Att 
% save(fnamesave,'res','s1','Wseed','t1','t2','Reps','testp','param','trial')
% 
