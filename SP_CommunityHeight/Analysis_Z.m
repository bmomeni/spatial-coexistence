%% Script to analyze the data from the cluster (spatial code)
% SKD 1/9/2020 code to analze output from cluster
% Compare community sizes between conditions


%% load/Rename Variables Z=.3 10/90
load ('CoexistenceCmp_WM_vs_SPYL_ExMT4_ABE_fp10_DMed0.018_TID1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_bt10_at15_Nc10_Nm5_qp50_qc50_Z0.3_Nz100_Nr10_Ns500_rndseed3725.mat')
DCell_AB5e8 = DCell;

NE0S_AB3 = NE0S;
NE0M_AB3 = NE0M;

V0MT_AB3 = V0MT;
V0ST_AB3 = V0ST;

DisAST_3 = DisAST;
DisBST_3 = DisBST;
DisEST_3 = DisEST;

ACmpA_3 = ACmpA;
ACmpB_3 = ACmpB;
ACmpE_3 = ACmpE;

CmpAST_3 = CmpAST;
CmpBST_3 = CmpBST;
CmpEST_3 = CmpEST;

CmpAMT_3 = CmpAMT;
CmpBMT_3 = CmpBMT;
CmpEMT_3 = CmpEMT;

%% load/Rename Variables Z=.5 10/90
load ('CoexistenceCmp_WMv2_vs_SPYL_ExMT4_ABE_fp10_DCell5e-08_TID1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_bt10_at15_Nc10_Nm5_qp50_qc50_Z0.5_Nz100_Nr10_Ns500_rndseed3725.mat')
DCell_FG5 = DCell;

NE0S_AB5 = NE0S;
NE0M_AB5 = NE0M;

V0MT_AB5 = V0MT;
V0ST_AB5 = V0ST;

DisBST_5 = DisBST;
DisEST_5 = DisEST;

ACmpB_5 = ACmpB;
ACmpE_5 = ACmpE;

CmpBST_5 = CmpBST;
CmpEST_5 = CmpEST;

CmpBMT_5 = CmpBMT;
CmpEMT_5 = CmpEMT;

%% load/Rename Variables Z=.7 10/90
load ('CoexistenceCmp_WM_vs_SPYL_ExMT4_ABE_fp10_DMed0.018_TID1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_bt10_at15_Nc10_Nm5_qp50_qc50_Z0.7_Nz100_Nr10_Ns500_rndseed3725.mat')
DCell_AB7 = DCell;

NE0S_AB7 = NE0S;
NE0M_AB7 = NE0M;

V0MT_AB7 = V0MT;
V0ST_AB7 = V0ST;

DisAST_7 = DisBST;
DisBST_7 = DisEST;

ACmpA_7 = ACmpB;
ACmpB_7 = ACmpE;

CmpAST_7 = CmpBST;
CmpBST_7 = CmpEST;

CmpAMT_7 = CmpBMT;
CmpBMT_7 = CmpEMT;

%% load/Rename Variables Z=.9 10/90
load ('CoexistenceCmp_WM_vs_SPYL_ExMT4_ABE_fp10_DMed0.018_TID1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_bt10_at15_Nc10_Nm5_qp50_qc50_Z0.3_Nz100_Nr10_Ns500_rndseed3725.mat')
DCell_AB9 = DCell;

NE0S_AB9 = NE0S;
NE0M_AB9 = NE0M;

V0MT_AB9 = V0MT;
V0ST_AB9 = V0ST;

DisAST_9 = DisBST;
DisBST_9 = DisEST;

ACmpA_9 = ACmpB;
ACmpB_9 = ACmpE;

CmpAST_9 = CmpBST;
CmpBST_9 = CmpEST;

CmpAMT_9 = CmpBMT;
CmpBMT_9 = CmpEMT;


%% Bootstrapping
%diffusion 5e8 Spatial
yD = NE0S_AB3(1,:);
MERDm_3A = mean(yD);
MERDci_3A = bootci(100, @mean, yD);
yD = NE0S_AB5(1,:);
MERDm_5A = mean(yD);
MERDci_5A = bootci(100, @mean, yD);
yD = NE0S_AB7(1,:);
MERDm_7A = mean(yD);
MERDci_7A = bootci(100, @mean, yD);
yD = NE0S_AB9(1,:);
MERDm_9A = mean(yD);
MERDci_9A = bootci(100, @mean, yD);

yD = NE0S_AB3(2,:);
MERDm_3B = mean(yD);
MERDci_3B = bootci(100, @mean, yD);
yD = NE0S_AB5(2,:);
MERDm_5B = mean(yD);
MERDci_5B = bootci(100, @mean, yD);
yD = NE0S_AB7(2,:);
MERDm_7B = mean(yD);
MERDci_7B = bootci(100, @mean, yD);
yD = NE0S_AB9(2,:);
MERDm_9B = mean(yD);
MERDci_9B = bootci(100, @mean, yD);

yD = NE0S_AB3(3,:);
MERDm_3E = mean(yD);
MERDci_3E = bootci(100, @mean, yD);
yD = NE0S_AB5(3,:);
MERDm_5E = mean(yD);
MERDci_5E = bootci(100, @mean, yD);
yD = NE0S_AB7(3,:);
MERDm_7E = mean(yD);
MERDci_7E = bootci(100, @mean, yD);
yD = NE0S_AB9(3,:);
MERDm_9E = mean(yD);
MERDci_9E = bootci(100, @mean, yD);

mean_Z_allA = [MERDm_3A,MERDm_5A,MERDm_7A,MERDm_9A];
range_Z_allA = [MERDci_3A,MERDci_5A,MERDci_7A,MERDci_9A];
mean_Z_allB = [MERDm_3B,MERDm_5B,MERDm_7B,MERDm_9B];
range_Z_allB = [MERDci_3B,MERDci_5B,MERDci_7B,MERDci_9B];
mean_Z_allE = [MERDm_3E,MERDm_5E,MERDm_7E,MERDm_9E];
range_Z_allE = [MERDci_3E,MERDci_5E,MERDci_7E,MERDci_9E];

%diffusion 5e-8
figure
hold on
errorbar([.3 .5 .7 .9], mean_Z_allA, mean_Z_allA - range_Z_allA(1,:), range_Z_allA(2,:) - mean_Z_allA, 's', 'Markersize', 8, 'color', [0 0.6 0]);
errorbar([.3 .5 .7 .9], mean_Z_allB, mean_Z_allB - range_Z_allB(1,:), range_Z_allB(2,:) - mean_Z_allB, 's', 'Markersize', 8, 'color', [1 0.4 1]);
errorbar([.3 .5 .7 .9], mean_Z_allE, mean_Z_allE - range_Z_allE(1,:), range_Z_allE(2,:) - mean_Z_allE, 's', 'Markersize', 8, 'color', [0 0 0]);
%ylim([0 3])
%xlim([0 12])
ylabel('Mean Richness')
xlabel('\it{Z} (cm)')

