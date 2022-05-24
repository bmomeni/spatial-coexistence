function [Ne, Cmp, Ne0, Cmp0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,cellRatio,rint,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau)

%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% 5: populations growth rates < 0.9 of max rate in the past 20 generations thrown out from coexistence
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

% rint: % matrix of interaction coefficients

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

%% Initial state 
% cellRatioArray = 1 / nCellType * ones(1,nCellType) % cell distrbution
cMed = zeros(Nm,1); % concentrations of interaction mediators

%% Cell-growth time-course
taurng = tau0:dtau:tauf;

nCell = TID * cellRatio'; % initial number of each cell type

for iRound = 1 : Nr
    cMed = TID / sum(nCell) * cMed;
    nCell = TID * cellRatio'; % initial number of each cell type
    nCell0 = nCell;
    %pause
    
    tau0 = 0; % in hours
    tau = tau0;
    
    count = 0;
    while (tau<=tauf-dtau) && (sum(nCell)<DilTh)
        
        count = count+1;
        tau = taurng(count);
        
        cMed = cMed + dtau*(B*nCell - A*nCell);
        cMed(cMed<0) = 0;
        
        reff = r0 + ((rint<0).*rint) * (cMed ./ kSat) + ((rint>=0).*rint) * (cMed ./ (cMed + kSat));
                    
        nCell = nCell + dtau * (reff .* nCell);
        nCell(nCell < ExtTh) = 0;
        
    end
    cellRatio = 1/sum(nCell)*nCell';
    
end
indx = 1:Nc;
Ne0 = indx(nCell>0);
Cmp0 = cellRatio(Ne0);
% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp0) > 0
    Cmp_sum = zeros(1,size(Cmp0,2));
    Cmp_sum(1,:) = sum(Cmp0);

    Cmp0 = Cmp0./Cmp_sum;
end

nGen = log(DilTh/TID)/log(2);
r = (nCell./nCell0).^(20/nGen);
stp = (r > abs(0.9*max(r)));
Ne = indx(stp);
Cmp = cellRatio(Ne);
% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp) > 0
    Cmp_sum = zeros(1,size(Cmp,2));
    Cmp_sum(1,:) = sum(Cmp);

    Cmp = Cmp./Cmp_sum;
end

return;
