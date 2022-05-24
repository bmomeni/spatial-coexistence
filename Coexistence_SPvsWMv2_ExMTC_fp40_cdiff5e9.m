%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
rndseed0 = 3725;
rng(rndseed0,'twister');

%number of ~interations
Ns = 500; % # of samples being screened
rndseed = round(100*Ns*rand(1,Ns));

Nc = 10; % # of cell types in the initial pool
Nm = 5; % # of mediators
TID = 1e4; % total initial cell density
kSat = 1e4; % interaction strength saturation level of each population %why saturation  -phrasing
ExtTh = 0.1; % population extinction threshold %high threshold to maintain (try lower); try number of individuals. (not part of pop)
DilTh = 1e7; % coculture dilution threshold

min_gr = 0.1; %min value for any growth rate (1/hr)
range_gr = 0.1; %range of variation of growth rate (0.1 - 0.2 currently)

ri0 = 0.2; % maximum interaction strength, 1/hr
fpi = 0.4; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
del = 0; % avg. decay rate (per hour)
at = 0.15; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.5; % probability of production link per population
qc = 0.5; % probability of influence link per population

%previously I had 200 = nGen
nGen = 100; % total number of generations of community growth simulated
GenPerRound = log(DilTh/TID)/log(2);
Nr = round(nGen/GenPerRound); % number of rounds of propagation
 
%% Diffusion parameters
DCell = 5e-9; % diffusion constant, cm^2/hour
DMed = 5*(10.^-5.5)*3600; % diffusion constant, cm^2/hour
% kY = 1e12; % total allowed yield % I had kY = 1e9*dz/5e-3; % total allowed yield (cells/dz)


%% Simulation domain
Z = 0.5; % community height in cm
Nz = 100;
dz = Z/(Nz-1);
dt = 0.1*(dz^2)/DMed;
dt = dtau/round(dtau/dt);
dc = 0.1*(dz^2)/DCell;
dc = round(dc/dt)*dt; %how often diffusion for cells takes place

kY = 1e9*dz/5e-3; % total allowed yield (cells/dz)

r0T = zeros(Nc,Ns); %matrix of zeros: 4 rows because 4 species and 1000 columns because 1000 samples being screened
SiT = zeros(Nm,Ns);

%Mediator matrix
%remains fixed for whole run
%varies between samples
AT = zeros(Nm,Nc,Ns); % Consumption
BT = zeros(Nm,Nc,Ns); % Production
DT = zeros(Nm,Ns); % Decay


%Cmp = Composition  M = mixed
%final composition
CmpAMT = zeros(Nc,Ns); %50/50
CmpBMT = zeros(Nc,Ns); %more Neg
CmpEMT = zeros(Nc,Ns); %More Pos


%Cmp = Composition  S = Spatial
%final composition
CmpAST = zeros(Nc,Ns);
CmpBST = zeros(Nc,Ns);
CmpEST = zeros(Nc,Ns);

%Interaction matrix
%fixed throughout
%different per Sample
rintAT = zeros(Nm,Nc,Ns); %50/50
rintBT = zeros(Nm,Nc,Ns); %More neg
rintET = zeros(Nm,Nc,Ns); %More Pos

%Viability
% Survival = 1
%Death/absense = 0
%at end
%M = Mixed
%S = Spatial
%ABE = rows
V0MT = zeros(3,Nc,Ns);
V0ST = zeros(3,Nc,Ns);

%Number of species that survive
%ABE = rows
%M = Mixed
%S = Spatial
NE0M = zeros(3,Ns);
NE0S = zeros(3,Ns);

% Distribution of cells
%need to name time?
dist = zeros(Nz,Nc,Nr);

DisAST = zeros(Nz,Nc,Nr,Ns);
DisBST = zeros(Nz,Nc,Nr,Ns);
DisEST = zeros(Nz,Nc,Nr,Ns);

%Composition of final community (no cut off values)
ACmpA = zeros(Nc, Ns);
ACmpB = zeros(Nc, Ns);
ACmpE = zeros(Nc, Ns);

for ns = 1:Ns
    disp(ns)
    tic
    
    rng(rndseed(ns),'twister');
    r0 = min_gr +range_gr*rand(Nc,1); % population reproduction rates, per hour
    %r0 = 0.1+0.1*rand(Nc,1); % population reproduction rates, per hour
    kSatVector = kSat * (0.5 + rand(Nm, 1)); % population levels for influence saturation
    
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
    
    rIntMatA = R .* DistInteractionStrengthMT_PA(Nc, Nm, ri0); % matrix of interaction coefficients, 50/50
    rIntMatB = R .* DistInteractionStrengthMT_PB(Nc, Nm, ri0, fpi); % matrix of interaction coefficients, more negative
    rIntMatE = R .* DistInteractionStrengthMT_PB(Nc, Nm, ri0, 1-fpi); % matrix of interaction coefficients, more positive
    
    % initial distribution for well-mixed communities
    PopDist = 1 / Nc * ones(1,Nc); % cell distribution; population ratios
    
    % initial distribution for spatial communities
    % SpPopDist = ones(Nz,Nc);
    % SpPopDist = Nz/sum(sum(SpPopDist))*SpPopDist;
    
    SpPopDist = zeros(Nz,Nc);
    for cdi = 1:Nc
        SpPopDist(floor((cdi-1)/Nc*Nz)+1:floor(cdi/Nc*Nz),cdi) = 1;
    end
    SpPopDist = Nz/sum(sum(SpPopDist))*SpPopDist;
    
    %% Simulating dynamics, M: Well-mixed
    [NeAM, CmpAM, NeAM0, CmpAM0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDist,rIntMatA,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
    [NeBM, CmpBM, NeBM0, CmpBM0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDist,rIntMatB,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
    [NeEM, CmpEM, NeEM0, CmpEM0] = WellmixedInteraction_DpMM_ExMTC_v2(Nr,r0,PopDist,rIntMatE,TID,A,B,kSat,ExtTh,DilTh,tauf,dtau);
    
    %% Simulating dynamics, M: Well-mixed
    [NeAS, CmpAS, NeAS0, CmpAS0, distA, AllCmpA] = Spatial1DInteraction_DpMM_ExMTC_SKD(Nr,r0,SpPopDist,rIntMatA,TID,A,B,kSat,kY,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc);
    [NeBS, CmpBS, NeBS0, CmpBS0, distB, AllCmpB] = Spatial1DInteraction_DpMM_ExMTC_SKD(Nr,r0,SpPopDist,rIntMatB,TID,A,B,kSat,kY,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc);
    [NeES, CmpES, NeES0, CmpES0, distE, AllCmpE] = Spatial1DInteraction_DpMM_ExMTC_SKD(Nr,r0,SpPopDist,rIntMatE,TID,A,B,kSat,kY,ExtTh/Nz,DilTh,tauf,dtau,Nz,Z,DCell,DMed,dt,dist,dc);
        
    V0AM = zeros(1,Nc);
    V0BM = zeros(1,Nc);
    V0EM = zeros(1,Nc);
    V0AM(NeAM) = 1;
    V0BM(NeBM) = 1;
    V0EM(NeEM) = 1;
    
    V0AS = zeros(1,Nc);
    V0BS = zeros(1,Nc);
    V0ES = zeros(1,Nc);
    V0AS(NeAS) = 1;
    V0BS(NeBS) = 1;
    V0ES(NeES) = 1;
    
    NE0M(:,ns) = sum([V0AM; V0BM; V0EM],2);
    NE0S(:,ns) = sum([V0AS; V0BS; V0ES],2);
    
    Cmp0AM = zeros(1,Nc);
    Cmp0BM = zeros(1,Nc);
    Cmp0EM = zeros(1,Nc);
    Cmp0AM(NeAM) = CmpAM;
    Cmp0BM(NeBM) = CmpBM;
    Cmp0EM(NeEM) = CmpEM;
    
    Cmp0AS = zeros(1,Nc);
    Cmp0BS = zeros(1,Nc);
    Cmp0ES = zeros(1,Nc);
    Cmp0AS(NeAS) = CmpAS;
    Cmp0BS(NeBS) = CmpBS;
    Cmp0ES(NeES) = CmpES;
    
    CmpAMT(:,ns) = Cmp0AM;
    CmpBMT(:,ns) = Cmp0BM;
    CmpEMT(:,ns) = Cmp0EM;
    CmpAST(:,ns) = Cmp0AS;
    CmpBST(:,ns) = Cmp0BS;
    CmpEST(:,ns) = Cmp0ES;
    
    r0T(:,ns) = r0;
    SiT(:,ns) = kSatVector;
    AT(:,:,ns) = A;
    BT(:,:,ns) = B;
    rintAT(:,:,ns) = rIntMatA';
    rintBT(:,:,ns) = rIntMatB';
    rintET(:,:,ns) = rIntMatE';
    V0MT(:,:,ns) = [V0AM; V0BM; V0EM];
    V0ST(:,:,ns) = [V0AS; V0BS; V0ES];
    
    DisAST (:,:,:,ns) = distA;
    DisBST (:,:,:,ns) = distB;
    DisEST (:,:,:,ns) = distE;
    
    ACmpA (:,ns) = AllCmpA;
    ACmpB (:,ns) = AllCmpB;
    ACmpE (:,ns) = AllCmpE;
    
    toc
end

save(strcat('CoexistenceCmp_WMv2_vs_SPYL_ExMT4_ABE_fp',num2str(round(100*fpi)),'_DMed',num2str(DMed,2),'_TID',num2str(TID,2),'_DilTh',num2str(DilTh,2),'_ExtTh',num2str(ExtTh),'_Ksat',num2str(kSat),'_ri',num2str(round(100*ri0)),'_bt',num2str(round(100*bt)),'_at',num2str(round(100*at)),'_Nc',num2str(Nc),'_Nm',num2str(Nm),'_qp',num2str(round(100*qp)),'_qc',num2str(round(100*qc)),'_Z',num2str(Z),'_Nz',num2str(Nz),'_Nr',num2str(Nr),'_Ns',num2str(Ns),'_rndseed',num2str(rndseed0),'.mat'))
