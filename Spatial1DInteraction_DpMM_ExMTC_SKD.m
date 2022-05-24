function [Ne, Cmp, Ne0, Cmp0, dist, AllCmp] = Spatial1DInteraction_DpMM_ExMTC_kY(Nr,r0,SpPopDist,rint,TID,A,B,kSat,kY,ExtTh,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc)

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

bb = 0;
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
%code I was using:  dz = Z/Nz;
dz = Z/(Nz-1);
%% Initial state
% CellDist = 1 / nCellType * ones(Nz,nCellType) % initial cell distrbution along z
cMed = zeros(Nz,Nm); % concentrations of interaction mediators at different heights
%% Cell-growth time-course
taurng = tau0:dtau:tauf;
trng = 0:dt:dtau;

nCell = TID * SpPopDist; % initial number of each cell type

for iRound = 1:Nr
    cMed = TID/sum(sum(1/Nz*nCell)) * cMed;
    nCell = TID * SpPopDist; % initial number of each cell type, Nz*Nc
    nCell0 = nCell; % initial number of each cell type
    
    tau0 = 0; % in hours
    tau = tau0;
    count = 0;
    
    while (tau<=tauf-dtau) && (sum(sum(nCell))<(Nz*DilTh))
        %I had: while (tau<=tauf-dtau) && (sum(sum(nCell))<DilTh)
        
        count = count+1;
        tau = taurng(count);
        
        % 1D diffusion finite difference matrix
        % Time lapse for diffusion
        for t = trng
            
            fMedm = dt*DMed/(dz^2)*[cMed(1,:); cMed(1:Nz-1,:)];
            fMedp = dt*DMed/(dz^2)*[cMed(2:Nz,:); cMed(Nz,:)];
            fMedc = -2*dt*DMed/(dz^2)*cMed;
            cMedCons = dt*nCell*A';
            cMed = cMed + dt*nCell*B' - cMedCons + (fMedm+fMedc+fMedp);
            cMed = cMed.*(cMed>0);
            
        end
        
        % 1D population diffusion
        nCell = nCell + dtau*DCell/(dz^2)*([nCell(1,:); nCell(1:Nz-1,:)] - 2*nCell + [nCell(2:Nz,:); nCell(Nz,:)]);
        reff = ones(Nz,1)*r0' + (1/kSat*cMed)*((rint<0).*rint)' + (cMed./(cMed + kSat))*((rint>=0).*rint)';
        
        %nCell = nCell + dtau*(reff.*nCell);
        nCell = nCell + dtau*(reff.*(ones(Nz,Nc)-1/kY*sum(nCell,2)*ones(1,Nc)).*nCell);
        
        nCell(nCell < ExtTh) = 0;
        
    end
    bb = 1+bb;
    dist(:,:,bb) = nCell;
    
    SpPopDist = Nz/sum(sum(nCell))*nCell;
    %What I had: SpPopDist = 1/sum(sum(nCell))*nCell;
    
end
indx = 1:Nc;
Ne0 = indx(sum(nCell,1)>0);
Cmp0 = sum(SpPopDist(:,Ne0),1);
% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp0) > 0
    Cmp_sum = zeros(1,size(Cmp0,2));
    Cmp_sum(1,:) = sum(Cmp0);
    
    Cmp0 = Cmp0./Cmp_sum;
end

nGen = log(DilTh/TID)/log(2);
r = (sum(nCell,1)./sum(nCell0,1)).^(20/nGen);
stp = (r > abs(0.9*max(r)));
Ne = indx(stp);
Cmp = sum(SpPopDist(:,Ne),1);
% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp) > 0
    Cmp_sum = zeros(1,size(Cmp,2));
    Cmp_sum(1,:) = sum(Cmp);
    
    Cmp = Cmp./Cmp_sum;
end

%Composition of final community (no cut off values)
indx = 1:Nc;
Ne0 = indx(sum(nCell,1)>0);

sum_nCell = sum(nCell);
holder = zeros(1,Nc);

for dd = 1:Nc
    if sum_nCell(1, dd) > 0
        holder(1, dd) = 1;
    end
end

AllCmp = holder.';

return;
