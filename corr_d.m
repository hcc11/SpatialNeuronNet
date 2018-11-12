function [C_d,COV_d,Cbar,COVbar,daxis,rate1,rate2,var1,var2]=corr_d(s1,s2,N1,N2,dim,Nc) 
% s1, s2: spike times of neuron Pop1 & Pop2, respectively, 2xN
% N1, N2: total # of neurons in Pop1 & Pop2, respectively
% Nc: 1x2 # of sampled neurons, e.g. Nc=[500, 500];
% only neuron in the center square [.25, .75]x[.25, .75] are sampled 
% Cbar=[C22, C12, C11], mean correlations 
% C_d: [20x3] correlation for each distance (daxis) 
% COV_d, COVbar: same as C_d, Cbar for covariance 
% rate1, rate2: mean rate (Hz) of the Nc sampled neurons 
% var1, var2: variance of the Nc sampled neurons 

T=floor(min([max(s1(1,:)) max(s2(1,:))]));

s1=s1(:,s1(1,:)<=T);
s2=s2(:,s2(1,:)<=T);

switch dim
    case '1D'
        I1=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
        I2=transpose(unique(s2(2,:)));
        
        I1=I1(I1<=N1);
        I2=I2(I2<=N2);
        
        Ic1=randsample(I1,Nc(1));
        Ic2=randsample(I2,Nc(2));

    case '2D'
        N11=round(sqrt(N1));
        N21=round(sqrt(N2));
        if size(s1,1)==3
        s1(2,:)=(s1(2,:)-1)*N11+s1(3,:); % x=ceil((I)/Nx1), y=mod(I-1,Nx1)+1
        end
        if size(s2,1)==3
        s2(2,:)=(s2(2,:)-1)*N21+s2(3,:);
        end
        I1=transpose(unique(s1(2,:)));
        I2=transpose(unique(s2(2,:)));
        
        I1=I1(I1<=N1);
        I2=I2(I2<=N2&I2>0);
        
        Ix10=(ceil(I1/N11))/N11;
        Iy10=(mod((I1-1),N11)+1)/N11;
        
        Ix20=(ceil(I2/N21))/N21;
        Iy20=(mod((I2-1),N21)+1)/N21;
        I1=I1(Ix10<0.75 & Ix10>0.25 & Iy10<0.75 & Iy10>0.25);
        I2=I2(Ix20<0.75 & Ix20>0.25 & Iy20<0.75 & Iy20>0.25);
        
        Ic1=randsample(I1,Nc(1));
        Ic2=randsample(I2,Nc(2));
end

% compute spike counts using sliding window 
Tw=200; % sliding window size 
Tburn=1000; 
time=0:1:T;

re1=zeros(Nc(1),length(time));
re2=zeros(Nc(2),length(time));

for mm=1:Nc(1)
    re1(mm,:)=hist(s1(1,Ic1(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1(mm)+1/4),time)*1e3;
end
for mm=1:Nc(2)
    re2(mm,:)=hist(s2(1,Ic2(mm)-1/4<s2(2,:) & s2(2,:)<=Ic2(mm)+1/4),time)*1e3;
end
re2_s=imfilter(re2(:,Tburn+1:end),ones(1,Tw)/Tw);re2_s=re2_s(:,Tw/2-1:end-Tw/2);
re1_s=imfilter(re1(:,Tburn+1:end),ones(1,Tw)/Tw);re1_s=re1_s(:,Tw/2-1:end-Tw/2);

ind1=mean(re1_s,2)>2;
ind2=mean(re2_s,2)>2;
re1_s=re1_s(ind1,:);
re2_s=re2_s(ind2,:);

switch dim 
    case '1D'
        Ix1=Ic1(ind1)/N1; Nc(1)=length(Ix1);
        Ix2=Ic2(ind2)/N2; Nc(2)=length(Ix2);
        D = pdist2([Ix2;Ix1],[Ix2;Ix1],'euclidean');
        D(D>0.5)=1-D(D>0.5); % periodic boundary condition
    case '2D'
        Ix1=(ceil((Ic1(ind1))/N11))/N11;
        Iy1=(mod((Ic1(ind1)-1),N11)+1)/N11;
        Nc(1)=length(Ix1);
        Ix2=(ceil((Ic2(ind2))/N21))/N21;
        Iy2=(mod((Ic2(ind2)-1),N21)+1)/N21;
        Nc(2)=length(Ix2);
        D = pdist2([[Ix2;Ix1],[Iy2;Iy1]],[[Ix2;Ix1],[Iy2;Iy1]],'euclidean');
end

COV=cov([re2_s; re1_s]');
Var=diag(COV);
var2=Var(1:Nc(2));var1=Var(Nc(2)+1:end);
rate1=mean(re1_s,2);
rate2=mean(re2_s,2);

R = COV./sqrt(Var*Var'); 
% Dx = pdist2([Ix2;Ix1],[Ix2;Ix1],'euclidean');
% Dx(Dx>0.5)=1-Dx(Dx>0.5); % periodic boundary condition
% Dy = pdist2([Iy2;Iy1],[Iy2;Iy1],'euclidean');
% Dy(Dy>0.5)=1-Dy(Dy>0.5); % periodic boundary condition
% D=sqrt(Dx.^2+Dy.^2);

dmax=0.5;
dd=0.025;
U=triu(ones(size(R)),1);
Utemp=zeros(size(U));
Utemp(1:Nc(2), 1:Nc(2))=U(1:Nc(2), 1:Nc(2));  % index for corr. within recurrent layer (E & I)
[Crr, COVrr, cbar_rr, COVbar_rr, daxis]=corr_dist(R,COV,D,Utemp, dd, dmax);
Utemp=zeros(size(U));
Utemp(1:Nc(2), Nc(2)+1:end)=U(1:Nc(2), Nc(2)+1:end);  % index for corr. btwn R & F
[Crf, COVrf, cbar_rf, COVbar_rf, daxis]=corr_dist(R,COV,D,Utemp, dd, dmax);
Utemp=zeros(size(U));
Utemp(Nc(2)+1:end, Nc(2)+1:end)=U(Nc(2)+1:end, Nc(2)+1:end);  % index for corr. btwn F & F
[Cff, COVff, cbar_ff, COVbar_ff, daxis]=corr_dist(R,COV,D,Utemp, dd, dmax);

C_d=[Crr,Crf,Cff];
COV_d=[COVrr,COVrf,COVff];
Cbar=[cbar_rr, cbar_rf, cbar_ff];
COVbar=[COVbar_rr, COVbar_rf, COVbar_ff];
end

function [c, cov_m,cbar, COVbar,daxis]=corr_dist(R,COV,D,U, dd, dmax)
% sort corr by distance according to index matrix U
% U is of same size as R & D
R=R(U==1); 
D=D(U==1);
COV=COV(U==1);
% dmax=0.5; 
% dd=0.025;
daxis=0:dd:dmax; 
[n, ind] = histc(D,daxis); 
n=n(1:end-1); % discard last bin (d>dmax) 
c=zeros(length(n),1);
cov_m=zeros(length(n),1);
for k=1:length(n)
    c(k)=mean(R(ind==k));
    cov_m(k)=mean(COV(ind==k));
end
daxis=daxis(1:end-1)+dd/2;  % center pts 
cbar=mean(R(D<=dmax)); % average corr.
COVbar=mean(COV(D<=dmax));
end 



