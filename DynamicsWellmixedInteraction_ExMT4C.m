function [tc, Xc, Cc] = DynamicsWellmixedInteraction_ExMT4C(Nr,r0,cellRatio,rint,TID,kSat,A,B,ExtTh,DilTh,tauf,dtau)

%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% 5: populations growth rates < 0.9 of max rate in the past 20 generations thrown out from coexistence
% rndseed = 1389;
% rand('twister',rndseed)

% nCellType = 15; % # of cell types
% nMediator = 6; % # of mediators
% nRound = 15; % number of rounds of propagation
% r0 = 0.08+0.04*rand(Nc,1); % population reproduction rates, per hour
% nInitialCell = 1e4; % total initial cells
% kSatLevel = 1e7; % interaction strength saturation level of each population
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

%% Initial state
% cellRatioArray = 1 / nCellType * ones(1,nCellType) % cell distrbution
cMed = zeros(Nm,1); % concentrations of interaction mediators

%% Cell-growth time-course
taurng = tau0:dtau:tauf;
Ntau = length(taurng);

tc = zeros(1,Nr);
Xc = zeros(Nc,Nr);
Cc = zeros(Nm,Nr);

nCell = TID * cellRatio'; % initial number of each cell type

cct = 0;
for iRound = 1:Nr
    cMed = TID / sum(nCell) * cMed;
    nCell = TID * cellRatio'; % initial number of each cell type
    %pause
    %disp(cellRatio)
    tau0 = 0; % in hours
    tau = tau0;
    count = 0;
    
    while (tau<=tauf-dtau) && (sum(nCell)<DilTh)
        count = count+1;
        tau = taurng(count);
        
        cMed = cMed + dtau*(B*nCell - A*nCell);
        cMed(cMed<0) = 0;
        
        reff = r0 + ((rint<0).*rint) * (cMed ./ kSat) + ((rint>=0).*rint) * (cMed ./ (cMed + kSat));
        
        nCell = nCell + dtau*(reff.*nCell);
        nCell(nCell < ExtTh) = 0;
        cct = cct+1;
        nCellt(:,cct) = nCell;
        cMedt(:,cct) = cMed;
    end
    cellRatio = 1/sum(nCell)*nCell';
    %disp(cellRatio)
    if iRound == 1
        tc(iRound) = taurng(count);
    else
        tc(iRound) = tc(iRound-1) + taurng(count);
    end
    Xc(:,iRound) = nCell;
    Cc(:,iRound) = cMed;
    
end
figure
plot(1:cct,nCellt','-')
set(gca,'YScale','log')
xlabel('Time-step')
ylabel('Populations')
 
% figure
% plot(1:cct,cMedt',':')
% set(gca,'YScale','log')
% xlabel('Time-step')
% ylabel('Chemical conc.')

return;
