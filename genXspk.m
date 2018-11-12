function sx=genXspk(p,Nx,T)
% generate spike trains of Layer 1
% p: param struct 
% Nx: # of neurons in Layer 1
% T: total simulation time 
switch p.stim_type
    case 'LocalCorr'
        cx=p.cx; 
        rX=p.rX;
        Nx1=sqrt(Nx);
        kc=round(Nx*cx);
        nspikes=round(T*1.1*rX/cx);
        tempspikes=sort(rand(nspikes,1)*(T));  % uniform distribution
        sx=zeros(2,kc*numel(tempspikes));
        for j=1:numel(tempspikes)
            sx(1,(j-1)*kc+1:j*kc)=tempspikes(j)+randn(kc,1)*p.taucorr;
            
            %%% Localized correlations
            spikeloc=rand(1,2);
            spikeindXs=mod(round((randn(kc,1)*p.sigmac+spikeloc(1))*Nx1)-1,Nx1);
            spikeindYs=mod(round((randn(kc,1)*p.sigmac+spikeloc(2))*Nx1)-1,Nx1);
            sx(2,(j-1)*kc+1:j*kc)=spikeindXs*Nx1+spikeindYs+1;
        end
        
    case 'Uncorr'
        % Generate uncorrelated spike trains for the feedforward layer
        sx=[];
        Nsource=p.Nsource; 
        rX=p.rX;
        for ns=1:Nsource
        tempspikes=cumsum(-log(rand(1,round(rX(ns)*(Nx/Nsource)*T*1.1)))/(rX(ns)*(Nx/Nsource)));
        tempspikes=tempspikes(tempspikes<T&tempspikes>0);
        sx_temp=zeros(2,numel(tempspikes));
        sx_temp(1,:)=tempspikes;
        sx_temp(2,:)=ceil(rand(1,size(sx_temp,2))*Nx/Nsource)+(ns-1)*Nx/Nsource;
        sx=[sx,sx_temp];
        end

    case 'GlobalCorr'
        sx=[];
        Nsource=p.Nsource; 
        rX=p.rX;cx=p.cx;
        for ns=1:Nsource
            % divid Nx for ns sources of correlation
            %%% globally correlated input
            kc_global=round(ceil(Nx/Nsource)*cx(ns));
            nspikes=round(T*1.1*rX(ns)/cx(ns));
%             tempspikes=sort(rand(nspikes,1)*(T-4*dt)+2*dt);  % uniform distribution
            tempspikes=cumsum(-log(rand(1,nspikes))/(rX(ns)/cx(ns)));
            sx_temp=zeros(2,kc_global*numel(tempspikes));
            for j=1:numel(tempspikes)
                sx_temp(1,(j-1)*kc_global+1:j*kc_global)=tempspikes(j)+randn(kc_global,1)*p.taucorr;
                
                %%% Global correlations
                spikeinds=randi(ceil(Nx/Nsource),1,kc_global)+ceil(Nx/Nsource)*(ns-1);
                sx_temp(2,(j-1)*kc_global+1:j*kc_global)=spikeinds;
            end
            sx=[sx, sx_temp];
        end
        
    case 'MultiSource' % when sigmaRx is not global 
        sx=[];
        Nsource=p.Nsource; 
        rX=p.rX;cx=p.cx;
        for ns=1:Nsource
            % divid Nx for ns sources of correlation
            %%% globally correlated input
            kc_global=round(ceil(Nx/Nsource)*cx(ns));
            nspikes=round(T*1.1*rX(ns)/cx(ns));
            tempspikes=cumsum(-log(rand(1,nspikes))/(rX(ns)/cx(ns)));
%             tempspikes=sort(rand(nspikes,1)*(T-4*dt)+2*dt);  % uniform distribution
            sx_temp=zeros(2,kc_global*numel(tempspikes));
            for j=1:numel(tempspikes)
                sx_temp(1,(j-1)*kc_global+1:j*kc_global)=tempspikes(j)+randn(kc_global,1)*p.taucorr;
                
                %%% Global correlations
                spikeinds=randi(ceil(Nx/Nsource),1,kc_global)*Nsource-(ns-1); 
                % mixed in space, mod(sx1,Nsource)=0,
                % mod(sx2,Nsource)=1, etc
                sx_temp(2,(j-1)*kc_global+1:j*kc_global)=spikeinds;
            end
            sx=[sx, sx_temp];
        end
    case 'spatialInput' % sigmac: spatial spread, centered at [.5 .5] 
         % peak rate is rX
%         fr=@(x,y) rX*exp(-((x-.5).^2+(y-.5).^2)/(2*sigmac^2))/(2*pi*sigmac^2);
%         %mean rate is rX
%         fr=@(x,y) rX*exp(-((x-.5).^2+(y-.5).^2));
%         center=[.5 .5]+[0.02, 0];
        rX=p.rX;
        center=p.center;
        Nx1=round(sqrt(Nx));
        CircRandN=@(mu,sigma,min,max,n)(mod(round(sigma*randn(n,1)+mu)-min,max-min+1)+min);
        mrate=rX*(2*pi*sigmac^2);
        tempspikes=cumsum(-log(rand(1,round(mrate*Nx*T*1.2)))/(mrate*Nx));
        tempspikes=tempspikes(tempspikes<T&tempspikes>0);
        sx=zeros(2,numel(tempspikes));
        sx(1,:)=tempspikes;
        Ix=CircRandN(center(1)*Nx1,(sigmac*Nx1),1,Nx1,numel(tempspikes));
        Iy=CircRandN(center(2)*Nx1,(sigmac*Nx1),1,Nx1,numel(tempspikes));
        sx(2,:)=(Ix-1)*Nx1+Iy;
end

[~,J]=sort(sx(1,:));
sx=sx(:,J);
sx=sx(:,sx(1,:)>0&sx(1,:)<=T);

