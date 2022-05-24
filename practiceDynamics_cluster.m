%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
simSize = 500;
NeAMs = zeros(1,simSize);
NeASs = zeros(1,simSize);
NeAMSPs = zeros(1,simSize);
NeBMs = zeros(1,simSize);
NeBSs = zeros(1,simSize);
NeBMSPs = zeros(1,simSize);
NeEMs = zeros(1,simSize);
NeESs = zeros(1,simSize);
NeEMSPs = zeros(1,simSize);



Nc = 10; % # of cell types in the initial pool
Nm = 5; % # of mediators
nInitialCell = 1e4; % total initial cells
TID = 1e4; % total initial cell density
kSat = 1e4; % interaction strength saturation level of each population
ExtTh = 0.1; % population extinction threshold
DilTh = 1e7; % coculture dilution threshold
ri0 = 0.2; % maximum interaction strength, 1/hr
fpi = 0.1; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
del = 0; % avg. decay rate (per hour)
at = 0.15; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.5; % probability of production link per population
qc = 0.5; % probability of influence link per population
kY = 1e9; % total allowed yield 1e12 is control %current results of 1/10/22 have 1e8

nGen = 100; % total number of generations of community growth simulated
rndseed0 = 3725;
rng(rndseed0,'twister');
rndseed = round(100*simSize*rand(1,simSize));
for i=1:simSize
disp(i)
rng(rndseed(i),'twister');
disp(rndseed);


GenPerRound = log(min(kY,DilTh)/TID)/log(2);
Nr = round(nGen/GenPerRound); % number of rounds of propagation

%% Diffusion parameters
DCell = 5e-9; % diffusion constant, cm^2/hour, equal to 5*10^-8 with constant added
DMed = 5*(10.^-6)*3600; % diffusion constant, cm^2/hour



%% Simulation domain
Z = 0.5; % community height in cm
Nz = 100;
dz = Z/(Nz-1);
dt = 0.1*(dz^2)/DMed;
%dt = dtau/round(dtau/dt);
dc = 0.1*(dz^2)/DCell;
dc = round(dc/dt)*dt; %how often diffusion for cells takes place


r0 = 0.1+0.1*rand(Nc,1); % population reproduction rates, per hour
kSatVector = kSat * (0.5 + rand(Nm, 1)); % population levels for influence saturation
dist = zeros(Nz,Nc,Nr);

%% Parameters
% Network configuration
% NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
% a probability q
R = NetworkConfig_Binomial(Nc,Nm,qc);
P = NetworkConfig_Binomial(Nc,Nm,qp);

% interaction matrix
alpha = at * (0.5+rand(Nc,Nm)); % consumption rates
beta = bt * (0.5+rand(Nc,Nm)); % mediator release rates
delta = del * (0.5+rand(Nm,1)); % decay rates, reusable chemicals only
A = (R.*alpha)';
B = (P.*beta)';

rIntA = R .* DistInteractionStrengthMT_PA(Nc, Nm, ri0); % matrix of interaction coefficients, 50/50
rIntB = R .* DistInteractionStrengthMT_PB(Nc, Nm, ri0, fpi); % matrix of interaction coefficients, more negative
rIntE = R .* DistInteractionStrengthMT_PB(Nc, Nm, ri0, 1-fpi); % matrix of interaction coefficients, more positive

%% Plot interaction network
rth = 0.001; % ignore interactions weaker than this threshold
cc = 11; % figure number for the plot
%PlotInteractionNetworkExMT4_StrengthBasalFitness(r0,rIntA,A,B,rth,cc)

% initial distribution for well-mixed communities
PopDist = 1 / Nc * ones(1,Nc); % cell distribution; population ratios

% initial distribution for spatial communities

CellDist = zeros(Nz,Nc);
for cdi = 1:Nc
    CellDist(floor((cdi-1)/Nc*Nz)+1:floor(cdi/Nc*Nz),cdi) = 1;
end
SpPopDist = 1/sum(sum(CellDist))*CellDist; %used to be 1

%[tcWM, XcWM, CcWM] = DynamicsWellmixedInteraction_ExMT4C(Nr,r0,PopDist,rIntA,TID,kSatVector,A,B,ExtTh,DilTh,tauf,dtau);
[NeAM, CmpAM, NeAM0, CmpAM0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDist,rIntA,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
[NeBM, CmpBM, NeBM0, CmpBM0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDist,rIntB,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
[NeEM, CmpEM, NeEM0, CmpEM0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDist,rIntE,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);

%[tcSP, XcSP, CcSP] = DynamicsSpatial1DInteraction_ExMTC(Nr,r0,SpPopDist,rIntA,TID,A,B,kSat,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt);
[NeAS, CmpAS, NeAS0, CmpAS0, distA, AllCmpA] = Spatial1DInteraction_DpMM_ExMTC_flexibleTimeStep(Nr,r0,SpPopDist,rIntA,TID,A,B,kSat,kY,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc);
[NeBS, CmpBS, NeBS0, CmpBS0, distB, AllCmpB] = Spatial1DInteraction_DpMM_ExMTC_flexibleTimeStep(Nr,r0,SpPopDist,rIntB,TID,A,B,kSat,kY,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc);
[NeES, CmpES, NeES0, CmpES0, distE, AllCmpE] = Spatial1DInteraction_DpMM_ExMTC_flexibleTimeStep(Nr,r0,SpPopDist,rIntE,TID,A,B,kSat,kY,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc);

PopDistSPA = zeros(size(PopDist));
PopDistSPB = zeros(size(PopDist));
PopDistSPE = zeros(size(PopDist));
PopDistSPA(NeAS) = CmpAS;
PopDistSPB(NeBS) = CmpBS;
PopDistSPE(NeES) = CmpES;

[NeAMSP, CmpAMSP, NeAM0SP, CmpAM0SP] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDistSPA,rIntA,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
[NeBMSP, CmpBMSP, NeBM0SP, CmpBM0SP] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDistSPB,rIntB,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
[NeEMSP, CmpEMSP, NeEM0SP, CmpEM0SP] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDistSPE,rIntE,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);

NeAMs(1,i) = length(NeAM);
NeASs(1,i) = length(NeAS);
NeBMs(1,i) = length(NeBM);
NeBSs(1,i) = length(NeBS);
NeEMs(1,i) = length(NeEM);
NeESs(1,i) = length(NeES);
NeAMSPs(1,i) = length(NeAMSP);
NeBMSPs(1,i) = length(NeBMSP);
NeEMSPs(1,i) = length(NeEMSP);
end

save(strcat('practiceDynamics_adjNr__simSize',num2str(simSize),'_dCell_5e-9_dMed5e-6*3600.mat'))

