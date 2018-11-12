function [Wrr,Wrf]=gen_weights(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,dimension)
% sigmaRR=[sigmaee, sigmaei; sigmaie, sigmaii]; 2x2 
% sigmaRX=[sigmaeX; sigmaiX]; 
% dimension='1D' or '2D'

% sort post syn index 

switch dimension
        
    case '2D'
        % x, y range: [1 Ne1], [1,Ni1], [1 Nx1] 
        % exc. ID [1, Ne], x=ceil(I/Ne1); y=(mod((I-1),Ne1)+1); I=(x-1)*Ne1+y 
        % inh. ID [Ne+1, Ne+Ni], x=ceil((I-Ne)/Ni1); y=(mod((I-Ne-1),Ni1)+1); I=(x-1)*Ni1+y+Ne; 
        
        % Connection widths
        sigmaeX=sigmaRX(1);
        sigmaiX=sigmaRX(2);
        sigmaee=sigmaRR(1,1);
        sigmaei=sigmaRR(1,2);
        sigmaie=sigmaRR(2,1);
        sigmaii=sigmaRR(2,2);
        
        pee0=Prr(1,1);
        pei0=Prr(1,2);
        pie0=Prr(2,1);
        pii0=Prr(2,2);
        pex0=Prx(1);
        pix0=Prx(2);
        
        Ne1=sqrt(Ne);
        Ni1=sqrt(Ni);
        Nx1=sqrt(Nx);
        betaee=sigmaee*(Ne1);
        betaei=sigmaei*(Ne1);
        betaie=sigmaie*(Ni1);
        betaii=sigmaii*(Ni1);
        betaex=sigmaeX*(Ne1);
        betaix=sigmaiX*(Ni1);
        % Kab is the total number of projections a neuron from
        % pop b makes to ALL neurons in pop a
        Kee=ceil(pee0*Ne);
        Kei=ceil(pei0*Ne);
        
        Kie=ceil(pie0*Ni);
        Kii=ceil(pii0*Ni);
        
        Kex=ceil(pex0*Ne);
        Kix=ceil(pix0*Ni);
        
        CircRandN=@(mu,sigma,min,max,n)(mod(round(sigma*randn(n,1)+mu)-min,max-min+1)+min);
        
        Ke=Kee+Kie; % Number of excitatory and inhibitory connections per cell
        Ki=Kei+Kii;
        Kx=Kex+Kix;
        Wrr=zeros(Ke*Ne+Ki*Ni,1,'int32'); % recurrent connections
        Wrf=zeros(Kx*Nx,1,'int32'); % feedforward connections
        
        for j=1:Ne
            % E pre, E post
            x_pre=ceil(j/Ne1);
            y_pre=mod(j-1,Ne1)+1;
            x_post=CircRandN(x_pre,betaee,1,Ne1,Kee);
            y_post=CircRandN(y_pre,betaee,1,Ne1,Kee);
            Wrr((1+(j-1)*Ke):(Kee+(j-1)*Ke))=sort((x_post-1)*Ne1+y_post);
            
            % E pre, I post
            x_pre=ceil(j/Ne1)*Ni1/Ne1;
            y_pre=(mod(j-1,Ne1)+1)*Ni1/Ne1;
            x_post=CircRandN(x_pre,betaie,1,Ni1,Kie);
            y_post=CircRandN(y_pre,betaie,1,Ni1,Kie);
            Wrr((Kee+1+(j-1)*Ke):(Kee+Kie+(j-1)*Ke))=sort((x_post-1)*Ni1+y_post+Ne);
        end
        for j=1:Ni
            % I pre, E post
            x_pre=ceil(j/Ni1)*Ne1/Ni1;
            y_pre=(mod(j-1,Ni1)+1)*Ne1/Ni1;
            x_post=CircRandN(x_pre,betaei,1,Ne1,Kei);
            y_post=CircRandN(y_pre,betaei,1,Ne1,Kei);
            Wrr((Ne*Ke+1+(j-1)*Ki):(Ne*Ke+Kei+(j-1)*Ki))=sort((x_post-1)*Ne1+y_post);
            % I pre, I post
            x_pre=ceil(j/Ni1);
            y_pre=(mod(j-1,Ni1)+1);
            x_post=CircRandN(x_pre,betaii,1,Ni1,Kii);
            y_post=CircRandN(y_pre,betaii,1,Ni1,Kii);
            Wrr((Ne*Ke+Kei+1+(j-1)*Ki):(Ne*Ke+Kei+Kii+(j-1)*Ki))=sort((x_post-1)*Ni1+y_post+Ne);
        end
        
        for j=1:Nx
            % X pre, E post
            x_pre=ceil(j/Nx1)*Ne1/Nx1;
            y_pre=(mod(j-1,Nx1)+1)*Ne1/Nx1;
            x_post=CircRandN(x_pre,betaex,1,Ne1,Kex);
            y_post=CircRandN(y_pre,betaex,1,Ne1,Kex);
            Wrf((Kx*(j-1)+1):(Kx*(j-1)+Kex))=sort((x_post-1)*Ne1+y_post);
            % X pre, I post
            x_pre=ceil(j/Nx1)*Ni1/Nx1;
            y_pre=(mod(j-1,Nx1)+1)*Ni1/Nx1;
            x_post=CircRandN(x_pre,betaix,1,Ni1,Kix);
            y_post=CircRandN(y_pre,betaix,1,Ni1,Kix);
            Wrf((Kex+1+(j-1)*Kx):(j*Kx))=sort((x_post-1)*Ni1+y_post+Ne);
        end
        
        case '1D'
        % Connection widths
        sigmaeX=sigmaRX(1);
        sigmaiX=sigmaRX(2);
        sigmaee=sigmaRR(1,1);
        sigmaei=sigmaRR(1,2);
        sigmaie=sigmaRR(2,1);
        sigmaii=sigmaRR(2,2);
        
        % Connection widths in units of neuron indices
        pee0=Prr(1,1);
        pei0=Prr(1,2);
        pie0=Prr(2,1);
        pii0=Prr(2,2);
        pex0=Prx(1);
        pix0=Prx(2);
        
        betaee=sigmaee*(Ne);
        betaei=sigmaei*(Ne);
        betaie=sigmaie*(Ni);
        betaii=sigmaii*(Ni);
        betaex=sigmaeX*(Ne);
        betaix=sigmaiX*(Ni);
        
        % Kab is the total number of projections a neuron from
        % pop b makes to ALL neurons in pop a
        Kee=ceil(pee0*Ne);
        Kei=ceil(pei0*Ne);
        
        Kie=ceil(pie0*Ni);
        Kii=ceil(pii0*Ni);
        
        Kex=ceil(pex0*Ne);
        Kix=ceil(pix0*Ni);
        
        CircRandN=@(mu,sigma,min,max,n)(mod(round(sigma*randn(n,1)+mu)-min,max-min+1)+min);
        
        Ke=Kee+Kie; % Number of excitatory and inhibitory connections per cell
        Ki=Kei+Kii;
        Kx=Kex+Kix;
        Wrr=zeros(Ke*Ne+Ki*Ni,1,'int32'); % recurrent connections
        Wrf=zeros(Kx*Nx,1,'int32'); % feedforward connections
        
        for j=1:Ne
            % E pre, E post
            Wrr((1+(j-1)*Ke):(Kee+(j-1)*Ke))=(CircRandN(j,betaee,1,Ne,Kee));
            
            % E pre, I post
            Wrr((Kee+1+(j-1)*Ke):(Kee+Kie+(j-1)*Ke))=(CircRandN(j/Ne*Ni+Ne,betaie,Ne+1,Ni+Ne,Kie));
        end
        for j=1:Ni
            % I pre, E post
            Wrr((Ne*Ke+1+(j-1)*Ki):(Ne*Ke+Kei+(j-1)*Ki))=(CircRandN(j/Ni*Ne,betaei,1,Ne,Kei));
            % I pre, I post
            Wrr((Ne*Ke+Kei+1+(j-1)*Ki):(Ne*Ke+Kei+Kii+(j-1)*Ki))=(CircRandN(j+Ne,betaii,Ne+1,Ne+Ni,Kii));
        end
        
        for j=1:Nx
            % X pre, E post
            Wrf((Kx*(j-1)+1):(Kx*(j-1)+Kex))=(CircRandN(j/Nx*Ne,betaex,1,Ne,Kex));
            % X pre, I post
            Wrf((Kex+1+(j-1)*Kx):(j*Kx))=(CircRandN(j/Nx*Ni+Ne,betaix,Ne+1,Ne+Ni,Kix));
        end
        
end
