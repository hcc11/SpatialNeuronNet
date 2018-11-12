% Simulations for factor analysis (Fig. 4E, Fig. 5B,C right, Fig. 6A,C and Fig. S6) 
% with fixed connectivity matrices 
% change the inhibitory time constant or projection width in Layer 3
% Inh='slow','fast' or 'broad'
% spike trains from Layer 2 are shared for the three inh conditions 
% run Inh='slow' first 

clear

data_folder='data/';

Wseed1_range=[8541; 45134; 48395; 3547; 14845;  71109; 99911; 98570;...
       68790; 16203 ];
   
Wseed2_range=[800281; 141887; 421762; 915736; 792208; 959493; ...
    157614;  970593; 957167; 485376]; 

testp.inE=[0];
testp.inI=[.2 .35];  % muI for Unatt. & att. conditions 

Inh='slow';  % Fig. 4E, Fig. 5B and Fig. 6A,C
% Inh='fast';  % Fig. S6A-E 
% Inh='broad';  % Fig. 5C, Fig. S6F-J 

%%%%%%%%%%%% for job array on cluster %%%%%%%%%%%%%%%% 
rng('shuffle');
AI = getenv('PBS_ARRAYID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
% job_dex range from 1 to Ntrial*Np*Nnws

Np=length(testp.inI)*length(testp.inE);  
Ntrial=15; 
Nnws=8; 

pid=mod(job_dex-1,Np)+1;
trial=mod(ceil(job_dex/Np)-1, Ntrial)+1;
nws=ceil(job_dex/(Np*Ntrial));  
ip1=1;
ip2=pid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

dim='2D';

opt.save=1; % save data 
opt.CompCorr=1; % compute correlations
    Nc=[500 500];  % # of E neurons sampled from Layer 3 & 2
if strcmp(Inh,'slow')
    opt.loadS1=0;  % generate spike trains from Layer 2
else
    opt.loadS1=1;  % use stored spike trains from Layer 2
end
opt.plotPopR=0; % plot population rate
opt.fixW=1; %  rng(Wseed1) before generating Wrr1, Wrf1; rng(Wseed2) before generating Wrr2, Wrf2
Wseed1=Wseed1_range(nws);
Wseed2=Wseed2_range(nws);

inE=testp.inE(ip1);
inI=testp.inI(ip2);
Jx=25*[1;0.6];
Prx=[ .05; .05];
dt=.01;
T=20000; % (ms)
s1_fname=strrep(sprintf('%sRF%s3layer_fixW_slowInh_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g_nws%.03g',...
    data_folder,dim,Jx(1),Jx(2),inE,inI,trial,dt,nws),'.','d'),

filename=strrep(sprintf('%sRF%s3layer_fixW_%sInh_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g_nws%.03g',...
    data_folder,dim,Inh,Jx(1),Jx(2),inE,inI,trial,dt,nws),'.','d'),

ParamChange={'param(2).Iapp', [inE;inI]; 'filename', filename;...
    'dt', dt; 'T', T; 'Nc',Nc;'param(2).Jx',Jx; 'param(2).Prx',Prx};

switch Inh 
    case 'slow' 
        ParamChange=[ParamChange;{'param(2).taudsyn(3,1)', 1;'param(2).taursyn(3,1)', 8}];
    case 'fast' 
        ParamChange=[ParamChange;{'param(2).taudsyn(3,1)', 1;'param(2).taursyn(3,1)', .5}];
    case 'broad'
        ParamChange=[ParamChange;{'param(2).sigmaRR', [.1 .2; .1 .2]}];
end

if opt.loadS1
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end
if opt.fixW
    ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
end

RF2D3layer(opt, ParamChange)

%% after all simulations, run spkcounts.m, then FA.m
