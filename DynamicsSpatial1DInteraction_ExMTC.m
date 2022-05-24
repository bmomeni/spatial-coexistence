function [tc, Xc, Cc] = DynamicsSpatial1DInteraction_ExMTC(Nr,r0,SpPopDist,rint,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt)

%% 1D spatial model for growth of interacting species
% Ex: explicitly including the mediators
% MT: multi-target mediators
% No concentration dependence for release or consumption
% rndseed = 1389;
% rand('twister',rndseed)

% Nc = 15; % # of cell types
% Nm = 6; % # of mediators
% Nr = 15; % number of rounds of propagation
% r0 = 0.08+0.04*rand(Nc,1); % population reproduction rates, per hour
% nInitialCell = 1e4; % total initial cells
% kSat = 1e7; % interaction strength saturation level of each population
% ExtTh = 0.1; % population extinction threshold
% DilTh = 1e10; % coculture dilution threshold
tau0 = 0;
% tauf = 250; % in hours
% dtau = 0.01; % in hours, cell growth update and uptake timescale
% at = 0.1; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
% bt = 1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
% mp = 3; % average number of production links per population
% mc = 2; % average number of consumption links per population

% intMat : % matrix of interaction coefficients

[Nc, Nm] = size(rint);

%% Parameters
% R = zeros(Nc,Nm);
% rndc = rand(Nc,Nm);
% R(rndc <= mc/Nm) = 1;
% P = zeros(Nc,Nm);
% rndp = rand(Nc,Nm);
% P(rndp <= mp/Nm) = 1;

% interaction matrix
% alpha = at*(0.5+rand(Nc,Nm)); % consumption rates
% beta = bt*(0.5+rand(Nc,Nm)); % mediator release rates
% A = (R.*alpha)';
% B = (P.*beta)';

%% Diffusion parameters
% DCell = 5e-8*3600; % diffusion constant, cm^2/hour
% DMed = 5e-6*3600; % diffusion constant, cm^2/hour

%% Simulation domain
% Z = 0.02; % community height in cm
% Nz = 50;
dz = Z/(Nz-1);


%% Initial state
% CellDist = 1 / nCellType * ones(Nz,nCellType) % initial cell distrbution along z
cMed = zeros(Nz,Nm); % concentrations of interaction mediators at different heights

%% Cell-growth time-course
taurng = tau0:dtau:tauf;
trng = 0:dt:dtau;

nCell = TID * SpPopDist;
SpPopDist = Nz * SpPopDist; % initial number of each cell type

tc = zeros(1,Nr);
Xc = zeros(Nc,Nr);
Cc = zeros(Nm,Nr);

cct = 0;
for iRound = 1:Nr
    cMed = TID/sum(sum(1/Nz*nCell)) * cMed;
    nCell = TID * SpPopDist; % initial number of each cell type, Nz*Nc
    %pause()
    
    tau0 = 0; % in hours
    tau = tau0;
    count = 0;
    
    while (tau<=tauf-dtau) && (sum(sum(nCell))<Nz*DilTh)
        count = count+1;
        tau = taurng(count);
        
        % 1D diffusion finite difference matrix
        % Time lapse for diffusion
        ct = 0;
        for t = trng
            ct = ct+1;
            
            fMedm = dt*DMed/(dz^2)*[cMed(1,:); cMed(1:Nz-1,:)];
            fMedp = dt*DMed/(dz^2)*[cMed(2:Nz,:); cMed(Nz,:)];
            fMedc = -2*dt*DMed/(dz^2)*cMed;
            cMedCons = dt*nCell*A';
            cMed = cMed + dt*nCell*B' - cMedCons + (fMedm+fMedc+fMedp);
            cMed = cMed.*(cMed>0);
        end
        
        % 1D population diffusion
        nCell = nCell + dtau*DCell/(dz^2)*([nCell(1,:); nCell(1:Nz-1,:)] - 2*nCell + [nCell(2:Nz,:); nCell(Nz,:)]);
        
        % growth rates per cell type per location, Nz*Nc
        reff = ones(Nz,1)*r0' + (1/kSat*cMed)*((rint<0).*rint)' + (cMed./(cMed + kSat))*((rint>=0).*rint)';
        
        nCell = nCell + dtau*reff.*nCell;
        % nCell = nCell + dtau*(reff.*(ones(Nz,Nc)-1/kY*sum(nCell,2)*ones(1,Nc)).*nCell);
        nCell(nCell < ExtTh) = 0;
        cct = cct+1;
        nCellst(:,cct) = 1/Nz*sum(nCell,1)';
        cMedst(:,cct) = 1/Nz*sum(cMed,1)';
        
    end
    SpPopDist = Nz/sum(sum(nCell))*nCell;
    %disp(SpPopDist)
    
    if iRound == 1
        tc(iRound) = taurng(count);
    else
        tc(iRound) = tc(iRound-1) + taurng(count);
    end
    Xc(:,iRound) = 1/Nz*sum(nCell,1);
    Cc(:,iRound) = 1/Nz*sum(cMed,1);
end
figure
plot(1:cct,nCellst')
set(gca,'YScale','log')
xlabel('Time-step')
ylabel('Populations')

figure
plot(1:cct,cMedst')
set(gca,'YScale','log')
xlabel('Time-step')
ylabel('Chemical conc.')

return;
