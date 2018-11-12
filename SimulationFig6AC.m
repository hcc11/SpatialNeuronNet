% run SimulationFig4E.m first 
% compute covariance components from Factor analysis as a function of distance 

addpath(genpath([pwd '/fa_Yu/']));
rng('shuffle');
dim='2D';
type='slowInh';

testp.inE=[0];
testp.inI=[.2 .35]; 
Jx=25*[1;0.6];
dt=0.01;

data_folder='data/';

Nsample=10; % number of sampling 
Np=2;
Ntrial=15;
M=30; % latent dimensionality to fit data 
Nc=500; % number of neuorns per sampling 
Ne1=200;
Ne=Ne1^2;  % number of neuorn in Layer 3 Exc. pop. 
numFolds=2; % number of cross-validation folds
Nws=8; % number of weight matrices realizations 
zDimList=M;
Tw=140;  % time window for spike counts 
Tburn=1000;
T=2e4; 
Nt=floor((T-Tburn)/Tw); % number of spike counts
T=Nt*Tw+Tburn;

fnamesave=strrep(sprintf('%sFA_model_%s_Jex%.03g_dist',data_folder,type,Jx(1)),'.','d'), 

%%%%%%%%%%% collect spike counts from Layer 3 Eec. pop. %%%%%%%%%%%%%
fname_sc=@(nws) strrep(sprintf('%sRF2D3layer_fixW_%s_Jex%.03g_Jix%.03g_inI_dt0d01_Nc%d_spkcount_nws%d',...
    data_folder,type, Jx(1),Jx(2),Ne,nws),'.','d'),
Ic2=1:Ne;
for nws=1:Nws
    for pid=1:Np
        E2{pid}=zeros(Ne, Nt*Ntrial);
        for trial=1:Ntrial
            inE=testp.inE(1);
            inI=testp.inI(pid); 
            filename=strrep(sprintf('%sRF%s3layer_fixW_%s_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f_dt%.03g_nws%.03g',...
                data_folder,dim,type,Jx(1),Jx(2),inE,inI,trial,dt,nws),'.','d'), % simulated with SimulationFig4E.m
            
            % compute spike counts using sliding window
            load(filename,'s2')
            s2=s2(:,s2(1,:)>=Tburn&s2(1,:)<=T&s2(2,:)<=Ne);
            s2(1,:)=s2(1,:)-Tburn;
            E2{pid}(:,(trial-1)*Nt+1:trial*Nt)=spktime2count(s2,Ic2,Tw,floor((T-Tburn)/Tw),1);        
        end
        
        if max(E2{pid}(:))<127
            E2{pid}=int8(E2{pid}); % spike counts are saved as int8 type
        else
            fprintf('max E2{%d} >127',pid) % exceed the range of int8 type
        end
        
    end
    save(fname_sc(nws),'E2','T','Tw','Tburn','filename','testp')
end


%%%%%%%%% compute Cov components vs. distance using Factor analysis %%%%%%%%%%%%%
for nws=1:Nws
    
    load(fname_sc(nws)) 
    FR_th=2; %firing rate threshold for sampling (Hz)
    Ic0=1:Ne;
    idx2=ones(Ne,1);
    for pid=1:Np
        idx2=idx2.*(mean(E2{pid},2)/Tw*1e3>FR_th);
    end
    for k=1:Nsample
        ss=(nws-1)*Nsample+k;
    Ic2=randsample(Ic0(idx2>0),Nc);

    Ix2=(ceil((Ic2)/Ne1))/Ne1; % x,y index for neural location 
    Iy2=(mod((Ic2-1),Ne1)+1)/Ne1;
    Dx = pdist2(Ix2',Ix2','euclidean');
    Dx(Dx>0.5)=1-Dx(Dx>0.5); % periodic boundary condition
    Dy = pdist2(Iy2',Iy2','euclidean');
    Dy(Dy>0.5)=1-Dy(Dy>0.5); % periodic boundary condition
    Dist=sqrt(Dx.^2+Dy.^2); % pairwise distance
    dmax=0.5;
    dd=0.025;
    daxis1=0:dd:0.7; % distance axis 
    U=triu(ones(size(Dist)),1);
    [n, ind] = histc(Dist(U==1),daxis1); % sort according to distance 
    n=n(1:end-1); % discard last bin (d>dmax)
    
    if nws==1&&k==1
        LL_d={zeros(length(n),M,Nsample*Nws),zeros(length(n),M,Nsample*Nws)};
        cov_d=zeros(length(n),Nsample*Nws,Np);
        LL_dxy={zeros(length(nx),length(ny),M,Nsample*Nws),zeros(length(nx),length(ny),M,Nsample*Nws)};
        Lambda=zeros(M,Nsample*Nws,Np);
    end
    
    for pid=1:Np 
        Y=double(E2{pid}(Ic2,:));
        COV=cov(Y');
        dim = crossvalidate_fa(Y, 'zDimList', zDimList,'showPlots',false,'numFolds',numFolds);
        
        L=dim(1).estParams.L;
        LL=L*L';
        [V,D]=eig(LL);
        la=diag(D);
        [la,I]=sort(la,'descend');
        Lambda(:,ss,pid)=la(1:M);
        
        for mm=1:M
            LL_tmp=V(:,I(mm))*V(:,I(mm))';
            LL_tmp=LL_tmp(U==1);  % covariance components
            
            for nn=1:length(n)
                LL_d{pid}(nn,mm,ss)=mean(LL_tmp(ind==nn));
            end
           
        end
        COV=COV(U==1);
        for nn=1:length(n)
            cov_d(nn,ss,pid)=mean(COV(ind==nn)); % total covariance 
        end
        
    end
    end
end
daxis1=daxis1(1:end-1)+dd/2;
save(fnamesave,'Tw','Lambda','LL_d','daxis1','cov_d','fname','M',...
    'Nc','numFolds','Nsample','testp')
