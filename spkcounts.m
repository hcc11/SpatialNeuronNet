%  run after Simulation_Fig4.m 
%   testp.inI=[.2 .35]; 
%   Inh='slow','fast' or 'broad';
%   filename=strrep(sprintf('%sRF2D3layer_fixW_%sInh_Jex25_Jix15_inE0_inI%.03g_ID%.0f_dt0d01_nws1',...
%         data_folder,Inh,inI,trial),'.','d');

rng('shuffle');

%%%%%%%%%%%% for job array on cluster %%%%%%%%%%%%%%%% 
AI = getenv('PBS_ARRAYID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
% job_dex range from 1 to Nnws
Nnws=8; 
nws=job_dex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Inh='slow'; 
% Inh='fast';  
% Inh='broad';   

dim='2D';
testp.inE=[0];
testp.inI=[.2 .35]; 
Jx=25*[1;0.6];
dt=0.01;

data_folder='data/';
 
Ntrial=15;
Np=length(testp.inI);
T=2e4;
Tw=140; %  window size
Tburn=1000;
Nt=floor((T-Tburn)/Tw); % # of trials
T=Nt*Tw+Tburn;
binedges=Tburn:Tw:T;
Ne1=200;
Nc=500; % sample size of neurons
FR_th=2; %firing rate threshold for sampling (Hz) 

fnamesave=sprintf('%sRF2D3layer_fixW_%sInh_Jex%.03g_Jix%.03g_inI_dt0d01_Nc%d_spkcount_nws%d',...
        data_folder,Inh, Jx(1),Jx(2),Nc,nws);
fnamesave=strrep(fnamesave,'.','d'),

%%%%%%%%% compute rates %%%%%%%%%%%%%
Ic2=1:Ne1^2;
for pid=1:Np
    rate=zeros(1,Ne1^2);
    for trial=1:Ntrial
        inE=testp.inE(1);
        inI=testp.inI(pid);
   
        filename=strrep(sprintf('%sRF%s3layer_fixW_%sInh_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g_nws%.03g',...
            data_folder,dim,Inh,Jx(1),Jx(2),inE,inI,trial,dt,nws),'.','d'),

        load(filename,'s2');
        s2=s2(:,s2(1,:)>Tburn&s2(1,:)<=T&s2(2,:)<=Ne1^2);
        rate=rate+hist(s2(2,:),Ic2)/(T-Tburn)*1e3;  
    end
    rate=rate/Ntrial;
    spkcount(pid).rate=rate; 
end
idx=ones(1,Ne1^2); 
for pid=1:Np
    idx=idx.*(spkcount(pid).rate>FR_th); 
end
Ic2=randsample(Ic2(idx>0),Nc);
clear idx rate; 
save(fnamesave,'spkcount','Ic2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(fnamesave,'Ic2')

Ic2=sort(Ic2); 
for pid=1:2
    spkcount(pid).Y=zeros(Nc, Nt*Ntrial);
    for trial=1:Ntrial
        inE=testp.inE(1);
        inI=testp.inI(pid);
    
        filename=strrep(sprintf('%sRF%s3layer_fixW_%sInh_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g_nws%.03g',...
            data_folder,dim,Inh,Jx(1),Jx(2),inE,inI,trial,dt,nws),'.','d'),
        load(filename,'s2');
        % compute spike counts using sliding window
        s2=s2(:,s2(1,:)>=Tburn&s2(1,:)<=T&s2(2,:)<=Ne1^2); 
        s2(1,:)=s2(1,:)-Tburn;
        re2=spktime2count(s2,Ic2,Tw,floor((T-Tburn)/Tw),0);
        spkcount(pid).Y(:,(trial-1)*Nt+1:trial*Nt)=re2; 
    end
end
save(fnamesave,'Tburn','T','Tw','spkcount','Ic2','testp')

