% codes for the neural field model in Fig 5A & B,C left 

%% Fig. 5A 
% provides wavenumber containing max real lam > 0 againt tau, sigma changes

clear 

data_folder='data/';
fnamesave=[data_folder 'neuralfield_stability'];

efit = @(mu)(mu.^2 .* (mu > 0));
ifit = @(mu)(mu.^2 .* (mu > 0));

Ne = 10000;
Ni = 10000;

Jee = 1;
Jei = -2;
Jie = 1.5;
Jii = -2.5;
pee = .008;
pei = .008;
pie = .008;
pii = .008;
KeeIn=pee*Ne;
KeiIn=(pei*Ni)*Ni/Ne;
KieIn=(pie*Ne)*Ne/Ni;
KiiIn=pii*Ni;

wee0=KeeIn*Jee;
wei0=KeiIn*Jei;
wie0=KieIn*Jie;
wii0=KiiIn*Jii;
W0 = [wee0 wei0;wie0 wii0];

Fin = [0.48  0.32];  % external inputs for [mu_e, mu_i]

reg = .012; % initial condition for firing rate re (kHz)
rig = .010; % initial condition for firing rate ri (kHz)

taue=5; % (ms)
taui=5;

sigmae=.1; % spatial connection width
sigmai=.1;

mu_i_range=[.32 .5]; % .32 for Unatt., .5 for Att. condition
for ss=1:2
    Fin(2) = mu_i_range(ss);
    % find steady state for rates
    q = .05;
    eps = 1e-6;
    change = 1;
    iter = 1;
    while (change > q*eps)
        ue = Fin(1) + wee0*reg + wei0*rig; % total current to e
        ui = Fin(2) + wie0*reg + wii0*rig; % total current to i
        
        renext = efit(ue);
        rinext = ifit(ui);
        changee = abs(renext-reg)/reg;
        changei = abs(rinext-rig)/rig;
        change = abs(changee + changei)/2;
        reg = (1-q)*reg + q*renext;
        rig = (1-q)*rig + q*rinext;
        iter = iter+1;
        if (iter > 50000)
            change = 0;
        end
    end
    
    % Calculate eigenvalues
    sigmai_use = .05:(.00125/2):.2;
    taui_use = 2.5:(.125/4):25;
    
    maxreallamplot=zeros(numel(sigmai_use),numel(taui_use));
    maxreallamindex=zeros(numel(sigmai_use),numel(taui_use));
    Fmodes=0:1:5;
    [f1,f2]=meshgrid(Fmodes,Fmodes);
    wavenums=sqrt(sort(unique(f1.^2+f2.^2)));
    
    for i = 1:numel(sigmai_use)
        for j = 1:numel(taui_use)
            maxreallam = zeros(1,numel(wavenums)); % eigenvalue w/ larger real part for each wave number
            % Fourier coefficients of W's, wn: wave number
            weet=@(wn)(wee0*exp(-2*wn.^2*pi^2*sigmae^2));
            weit=@(wn)(wei0*exp(-2*wn.^2*pi^2*sigmai_use(i)^2));
            wiet=@(wn)(wie0*exp(-2*wn.^2*pi^2*sigmae^2));
            wiit=@(wn)(wii0*exp(-2*wn.^2*pi^2*sigmai_use(i)^2));
            Wt=@(wn)([weet(wn) weit(wn); wiet(wn) wiit(wn)]);
            
            Tau_e = taue;
            Tau_i = taui_use(j);
            
            Tau = [-1/Tau_e 0; 0 -1/Tau_i];
            
            ge = 2*ue;
            gi = 2*ui;
            
            GTau=[ge/Tau_e 0; 0 gi/Tau_i];
            
            for k=1:numel(wavenums)
                Matrix = Tau + GTau*Wt(wavenums(k));
                maxreallam(k)=max(real(eig(Matrix)));
            end
            [maxreallamplot(i,j),maxreallamindex(i,j)] = max(real(maxreallam));
        end
    end
    
    maxReWavNum = maxreallamindex.*(maxreallamplot>=0); % 0 for negative eigenvalues
    maxReWavNum(maxReWavNum==0) = maxReWavNum(maxReWavNum==0) -1;  % convert from index to wn, -1 for negative eigenvalues
    maxReWavNum(maxReWavNum>0) = wavenums(maxReWavNum(maxReWavNum>0)); % convert from index to wn
    res(ss).maxReWavNum=maxReWavNum;
    res(ss).muI=Fin(2);

end

save(fnamesave,'res','sigmai_use','taui_use','wavenums','taue','Fin') 

%% Fig. 5B,C left
% eigenval (max real) vs wn, for various taui

fnamesave=[data_folder 'neuralfield_Lambda']; 

Fin = [0.48  0.32]; 

reg = .012;
rig = .010;

taue=5;
taui_use = 5:2.5:15;
sigmae=.1;

ss=0;
for sigmai=[.1 .15]
    ss=ss+1;
    weet=@(wn)(wee0*exp(-2*wn.^2*pi^2*sigmae^2));
    weit=@(wn)(wei0*exp(-2*wn.^2*pi^2*sigmai^2));
    wiet=@(wn)(wie0*exp(-2*wn.^2*pi^2*sigmae^2));
    wiit=@(wn)(wii0*exp(-2*wn.^2*pi^2*sigmai^2));
    Wt=@(wn)([weet(wn) weit(wn); wiet(wn) wiit(wn)]);
    
    % find steady state for rates
    q = .05;
    eps = 1e-6;
    change = 1;
    iter = 1;
    while (change > q*eps)
        ue = Fin(1) + wee0*reg + wei0*rig;
        ui = Fin(2) + wie0*reg + wii0*rig;
        
        renext = efit(ue);
        rinext = ifit(ui);
        changee = abs(renext-reg)/reg;
        changei = abs(rinext-rig)/rig;
        change = abs(changee + changei)/2;
        reg = (1-q)*reg + q*renext;
        rig = (1-q)*rig + q*rinext;
        iter = iter+1;
        if (iter > 50000)
            change = 0;
        end
    end
    % Calculate eigenvalues
    Fmodes=0:1:5;
    [f1,f2]=meshgrid(Fmodes,Fmodes);
    wavenums=sqrt(sort(unique(f1.^2+f2.^2)));
    maxreallam=zeros(numel(taui_use),numel(wavenums));
    
    for j = 1:1:numel(taui_use)        
        Tau_e = taue;
        Tau_i = taui_use(j);
        
        Tau = [-1/Tau_e 0; 0 -1/Tau_i];
        
        ge = 2*ue;
        gi = 2*ui;
        
        GTau=[ge/Tau_e 0; 0 gi/Tau_i];
        
        for k=1:numel(wavenums)
            Matrix = Tau + GTau*Wt(wavenums(k));
            maxreallam(j,k)=max(real(eig(Matrix)));
        end % k loop
    end % j loop
    Lambdas(ss).sigmai=sigmai;
    Lambdas(ss).maxreallam=maxreallam;
end

save(fnamesave,'Lambdas','wavenums','taui_use','taue','sigmae','Fin') 
