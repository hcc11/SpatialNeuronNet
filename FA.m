%  run after SimulationFig4E.m & spkcounts.m
% spike counts data provided for nws=1 (set Nnws=1) 

addpath(genpath([pwd '/fa_Yu/']));
clear

data_folder='data/';

Inh='slow'; fnamesave=[data_folder 'FA_model_Jex25.mat']; % Fig. 4b, Fig. 5b bottom (Nc=500, Nsample=1)
% Inh='fast'; fnamesave=[data_folder 'FA_model_fastInh_Jex25.mat']; % Fig. 4c 
% Inh='broad'; fnamesave=[data_folder 'FA_model__broadInh_Jex25.mat']; % Fig. 4d,  Fig. 5c bottom (Nc=500, Nsample=1)

dim='2D';

testp.inE=[0];
testp.inI=[.2 .35];
Jx=25*[1;0.6];
dt=0.01;

Nws=8; % number of weight matrix realizations
Nsample=80;  % 10 samples per weight matrix realization
Np=2;
M=5; % latent dimensionality to fit data 
Nc=50; % neuron number 
numFolds=2; % number of cross-validation folds
corr=zeros(Nsample,Np);  %mean correlation 
sumLL=zeros(M,Nsample,Np); % cross-validated log likelihood
LLopt=zeros(Nsample,Np); % mode that maximizes the cross-validated log likelihood
Lambda=zeros(M,Nsample,Np);  % eigenvalues 
COVm=zeros(Nsample,Np);  % mean covariance  (raw)
Qm=zeros(Nsample,Np);   % residual covarince
eigvector1=NaN(Nc,Nsample,Np);  % first eigenvector 
SIGN=zeros(Nsample,Np);  % sign of the mean of the first eigenvector 
zDimList=1:M;

for nws=1:Nws
    fname=sprintf('%sRF2D3layer_fixW_%sInh_Jex25_Jix15_inI_dt0d01_Nc500_spkcount_nws%d',data_folder,Inh,nws);
    data=load(fname);
    idx=randperm(Nc*10);
    
    for pid=1:2
        for k=1:10
            ss=(nws-1)*10+k;
            Y=data.spkcount(pid).Y(idx((k-1)*Nc+1:k*Nc),:);      
            COV=cov(Y');
            U=triu(ones(size(COV)),1);
            R=corrcov(COV);
            corr(ss,pid)=mean(R(U==1));
            dim = crossvalidate_fa(Y, 'zDimList', zDimList,'showPlots',false,'numFolds',numFolds);
            
            sumLL(:,ss,pid)=[dim.sumLL];
            LLopt(ss,pid) = find(sumLL(:,ss,pid) == max(sumLL(:,ss,pid)),1);
            
            L=dim(M).estParams.L;
            LL=L*L';
            [V,D]=eig(LL);
            la=diag(D);
            [m,I]=max(la);
            la=sort(la,'descend');
            Lambda(:,ss,pid)=la(1:M);
            eigvector1(:,ss,pid)=V(:,I)*sign(sum(V(:,I)));
            SIGN(ss,pid)=sign(sum(V(:,I)));
            
            L1=dim(1).estParams.L;
            LL1=L1*L1';
            Q=COV-LL1;
            COVm(ss,pid)=mean(COV(U==1));
            Qm(ss,pid)=mean(Q(U==1));
            
        end
    end
end
figure
colorAU= [0    0.5000    0.4000;
     0.9290    0.6940    0.1250];
state={'Unattended','Attended'};
hold on
for pid=1:Np
errorbar(1:M, mean(Lambda(:,:,pid),2),std(Lambda(:,:,pid),[],2)/sqrt(Nsample),'-','color',colorAU(pid,:))
text(.8, 1-pid*.1,state{pid},'unit','n','Horiz','center','color',colorAU(pid,:))
end
xlabel('eigenmode')
ylabel('eigenvalue')

save(fnamesave,'corr','COVm','Qm','LLopt','Lambda','eigvector1','fname','M','Nc','numFolds','Nsample','testp')


