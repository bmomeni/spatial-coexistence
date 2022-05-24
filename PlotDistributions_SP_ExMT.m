%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
rndseed0 = 3752;
%rand('twister',rndseed0);
rng(rndseed0,'twister');

Nc = 5; % # of cell types in the initial pool
Nm = 3; % # of mediators
%Nr = 10; % number of rounds
nInitialCell = 1e4; % total initial cells
TID = 1e4; % total insitial cell density
kSat = 1e4; % interaction strength saturation level of each population
ri0 = 0.2; % maximum interaction strength, 1/hr
fpi = 0.1; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
at = 0.15; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.5; % probability of production link per population
qc = 0.5; % probability of influence link per population
Z = 0.5; % Thickness of community in cm
Nz = 100; % discretization of height
kY = 1e9; % total allowed yield
DilTh = 1e7; % coculture dilution threshold
ExtTh = 0.1; % population extinction threshold

nGen = 100; % total number of generations of community growth simulated
GenPerRound = log(DilTh/TID)/log(2);
Nr = round(nGen/GenPerRound); % number of rounds of propagation

%% Diffusion parameters
DCell = 5*(10.^-9); % diffusion constant, cm^2/hour
DMed = 5*(10.^-6)*3600; % diffusion constant, cm^2/hour

%% Simulation domain
dz = Z/(Nz-1);
dt = 0.1*(dz^2)/DMed;
dt = dtau/round(dtau/dt);

r0 = 0.1+0.1*rand(Nc,1); % population reproduction rates, per hour
kSatVector = kSat * (0.5 + rand(Nm, 1)); % population levels for influence saturation
%% Parameters
% Network configuration
% NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
% a probability q
R = NetworkConfig_Binomial(Nc,Nm,qc);
P = NetworkConfig_Binomial(Nc,Nm,qp);

% consumption and release matrices
alpha = at * (0.5+rand(Nc,Nm)); % consumption rates
beta = bt * (0.5+rand(Nc,Nm)); % mediator release rates
A = (R.*alpha)'; % Nm*Nc
B = (P.*beta)'; % Nm*Nc

% interaction matrix
rIntA = R .* DistInteractionStrengthMT_PA(Nc, Nm, ri0); % matrix of interaction coefficients, 50/50
rIntB = R .* DistInteractionStrengthMT_PB(Nc, Nm, ri0, fpi); % matrix of interaction coefficients, more negative
rIntE = R .* DistInteractionStrengthMT_PB(Nc, Nm, ri0, 1-fpi); % matrix of interaction coefficients, 50/50

%% Initial state
% CellDist = 1 / nCellType * ones(Nz,nCellType) % initial cell distrbution along z
%cMed = zeros(Nz,Nm); % concentrations of interaction mediators at different heights

%% Cell-growth time-course
cMed = zeros(Nz,Nm); % concentrations of interaction mediators at different heights
taurng = tau0:dtau:tauf;
Nt = length(taurng);
trng = 0:dt:dtau;

% tc = zeros(1,Nr*Nt);
% Xc = zeros(Nz,Nc,Nr*Nt);
% Cc = zeros(Nz,Nm,Nr*Nt);

% spatial distributions are initially random
% CellDist = rand(Nz,Nc);
% CellDist = 1/sum(sum(CellDist))*CellDist;
CellDist = zeros(Nz,Nc);
for cdi = 1:Nc
    CellDist(floor((cdi-1)/Nc*Nz)+1:floor(cdi/Nc*Nz),cdi) = 1;
end
CellDist = 1/sum(sum(CellDist))*CellDist;

cct = 0;
%nCell = nInitialCell * CellDist; % initial number of each cell type
nCell = TID * Nz * CellDist; % initial number of each cell type, Nz*Nc
nCell0 = nCell; % initial number of each cell type

nCellt = zeros(Nc,Nt);
cMedt = zeros(Nm,Nt);
frames = struct('cdata',[],'colormap',[]);

tic
count = 0;
crt = 0;
crt = crt+1;
nCellt(:,crt) = sum(nCell,1);
cMedt(:,crt) = sum(cMed,1);
f = figure(2);
plot(1:Nz,nCell)
legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10')
xlabel('height')
%add good title for this graph and y-axis
width = 560;
height = 415;
set(f,'Position',[15 15 width height])
frames(crt) = getframe(f);
figure(3)
plot(1:Nz,cMed)
legend('C1','C2','C3','C4','C5')
xlabel('height')
    %pause(1)
for iRound = 1:Nr
    %cMed = nInitialCell/sum(sum(1/Nz*nCell)) * cMed;
    %nCell = nInitialCell * CellDist; % initial number of each cell type, Nz*Nc
    
    cMed = TID/sum(sum(1/Nz*nCell)) * cMed;
    nCell = TID * CellDist; % initial number of each cell type, Nz*Nc
    nCell0 = nCell;
    
    tau0 = 0; % in hours
    tau = tau0;
    %count = 0;
    
    while (tau<=tauf-dtau) && (sum(sum(nCell))<DilTh)%*1/Nz)
        
        count = count+1;
        tau = tau+dtau;
        %tau = taurng(count);
        
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
        reff = ones(Nz,1)*r0' + (1/kSat*cMed)*((rIntE<0).*rIntE)' + (cMed./(cMed + kSat))*((rIntE>=0).*rIntE)';
        
        %nCell = nCell + dtau*reff.*nCell;
        nCell = nCell + dtau*(reff.*(ones(Nz,Nc)-1/kY*sum(nCell,2)*ones(1,Nc)).*nCell);
        nCell(nCell < ExtTh) = 0;
        
%         if mod(count,1000)==1
%             crt = crt+1;
%             nCellt(:,crt) = sum(nCell,1);
%             cMedt(:,crt) = sum(cMed,1);
%             f = figure(2);
%             plot(1:Nz,nCell)
%             legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10')
%             xlabel('height')
%             %add good title for this graph and y-axis
%             width = 560;
%             height = 415;
%             set(f,'Position',[15 15 width height])
%             frames(crt) = getframe(f);
%             figure(3)
%             plot(1:Nz,cMed)
%             legend('C1','C2','C3','C4','C5')
%             xlabel('height')
%             %pause(1)
%         end
    end
    crt = crt+1;
    nCellt(:,crt) = sum(nCell,1);
    cMedt(:,crt) = sum(cMed,1);
    f = figure(2);
    plot(1:Nz,nCell)
    legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10')
    xlabel('height')
    %add good title for this graph and y-axis
    width = 560;
    height = 415;
    set(f,'Position',[15 15 width height])
    frames(crt) = getframe(f);
    figure(3)
    plot(1:Nz,cMed)
    legend('C1','C2','C3','C4','C5')
    xlabel('height')
    %pause(1)
    CellDist = Nz/sum(sum(nCell))*nCell;  %could be Nz as well, need to figure out this line
    
end
toc
nGen = log(DilTh/TID)/log(2);                                                    
r = (sum(nCell,1)./sum(nCell0,1)).^(20/nGen);
%r = sum(nCell,1)./sum(nCell0,1);
stp = (r > abs(0.9*max(r)));
indx = 1:Nc;
Ne = indx(stp);
Cmp = CellDist(:,Ne);

% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp) > 0
    Cmp_sum = zeros(1,size(Cmp,2));
    Cmp_sum(1,:) = sum(Cmp);
    
    Cmp = Cmp./Cmp_sum;
end

figure
semilogy(0:crt-1,nCellt(:,1:crt)')
xlabel('Number of Rounds')
ylabel('Population Size per Species')

figure
semilogy(0:crt-1,cMedt(:,1:crt)')
xlabel('Number of Rounds')
ylabel('Mediator Concentration')

vw = VideoWriter('spatial_test.m4v','MPEG-4');
vw.FrameRate = 3;
open(vw)
writeVideo(vw,frames)
close(vw)