%  demo code for two-layer network simulations (related to Fig. 3)
%  calls RF2D3layer.m for network simulation 

clear 

Wtype='broadRec';  taudsyni=8; % (ms)   % Fig. 3Aiv
% Wtype='broadRec';  taudsyni=0.5; % (ms)   % Fig. 3Aii
% Wtype='uniformW'; taudsyni=8; % Fig. 3Aiii
% Wtype='uniformW'; taudsyni=0.5; % Fig. 3Ai
data_folder='data/';   % folder name to store data 
dim='2D'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic


%%%%%%%% set options and parameters to change %%%%%%%%%%%
opt.save=1; % save data 
opt.CompCorr=1; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.Layer1only=1; % 1 for two-layer network, 0 for three-layer network  
opt.loadS1=0;
opt.plotPopR=1; % plot population rate
opt.fixW=0;  % use the same weight matrices for multiple simulations 
%     Wseed1=Wseed1_range(nws); % ransom seed for generating weight matrices
%     Wseed2=Wseed2_range(nws); 

dt=.05; % time step size for integration 
T=2000; % total simulation time (ms) 
filename=strrep(sprintf('%sRF%s2layer_%s_tausyni%.03g_test',...
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

%%%%%% plot correlation vs distance %%%%%%%%%%%%%
load(filename,'C','daxis')
figure
plot(daxis,C,'linewidth',1) 
hold on
colororder=lines(3); 
text(.7, 0.9,'Corr (E2-E2)','unit','n','Horiz','l','color',colororder(1,:))
text(.7, 0.8,'Corr (E2-E1)','unit','n','Horiz','l','color',colororder(2,:))
text(.7, 0.7,'Corr (E1-E1)','unit','n','Horiz','l','color',colororder(3,:))
xlabel('distance (a.u.)')
ylabel('correlation')

%%%%%%% raster movie %%%%%%%%%%%%%%%%
load(filename,'s1')
t1=0; % start time (ms)
t2=1000; % end time (ms)
raster2D_ani(s1,t1,t2,200)

