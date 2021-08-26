clear all;
close all;
clc;

%% Mesh convergence
% Euler-Bernoulli
Euler = 0.001190476190476;

% MATLAB

CST_dofs = [30,594,2222,3402,5022,8282,10302,12322,14762,17182,19602,22842,24462,25662];
CST = [1.489594775595639e-04,0.001064951716968,0.001190388895012,0.001200283991053,...
    0.001203717354782,0.001213720937214,0.001214289517394,0.001214598443731,0.001219601547756,0.001219790119117,...
    0.001219912491413,0.001223069785971,0.001224238175844,0.001223154500608];

LST_dofs = [90,450,1066,4050,15778,20930,26082];
LST = [0.001182631300959,0.001226797307065,0.001228936101706,0.001230944439497,0.001232875908865,0.001232897812417,0.001232906522802];


TR_dofs = [81,567,3843,8019,14883,18513,24633];
TR = [3.720323776394732e-04,9.187284504480370e-04,0.001184102288880,0.001205424706547,0.001220006630046,0.001220504038096,0.001225785463890];


% ABAQUS

LST_dofs_A = [90,594,2210,8514,13202,20930,26082];
LST_A = [0.00118263,0.00122808,0.00123029,0.00123223,0.00123285,0.0012329,0.00123291];

CST_dofs_A = [30,594,2210,5346,10302,20930];
CST_A = [0.000148959,0.00106495,0.00118384,0.00120405,0.00121429,0.00122499];


% Figures

% figure;
% semilogx(1./CST_dofs,CST./Euler,'*-',1./CST_dofs_A,CST_A./Euler,'o-',1./LST_dofs,...
%     LST./Euler,'*-',1./LST_dofs_A,LST_A./Euler,'o-',1./TR_dofs,TR./Euler,'*-');
% grid on;
% xlabel('1/DOFs'); ylabel('FE/EB');
% title('Midspan vertical displacement');
% legend('MATLAB CST','ABAQUS CST','MATLAB LST','ABAQUS LST','MATLAB TR');


%% Code-performance
truss_K = [0.024577,0.028822,0.031899,0.044881,0.085547,0.22814,0.90918,3.7748];
truss_inv = [0.00033109,0.00061051,0.0018556,0.0060828,0.024911,0.10484,0.52585,2.7079];
truss_optimisation = [0.034978,0.036621,0.038603,0.041353,0.047507,0.075832,0.23268,0.85006];

truss_K_MATLAB = [0.012772,0.012961,0.014516,0.036542,0.098375];
truss_inv_MATLAB = [0.0006354,0.0023271,0.031123,0.34077,8.5015];

CST_K = [0.025079,0.026446,0.029257,0.050515,0.10068,0.2546,1.0562,4.3942];
CST_inv = [0.00029118,0.0004006,0.0012842,0.0064867,0.031564,0.11104,0.52212,2.6641];
CST_optimisation = [0.038434,0.039004,0.041984,0.043423,0.045402,0.072949,0.26726,1.0051];

CST_K_MATLAB = [0.038388,0.040148,0.043277,0.082537,0.17222];
CST_inv_MATLAB = [0.00061869,0.0014643,0.01764,0.34659,8.6908];

LST_K = [0.026093,0.030714,0.047502,0.12869,0.2987,1.2434,5.3363,20.2126];
LST_inv = [0.00054692,0.0018272,0.0063568,0.074497,0.2974,3.466,17.0554,84.451];
LST_optimisation = [0.032132,0.033381,0.037232,0.043121,0.052531,0.1254,0.44629,1.9171];

LST_K_MATLAB = [0.53574,0.53807,0.63554,0.78885];
LST_inv_MATLAB = [0.0052043,0.012145,0.26813,9.9065];

NX = [20,40,80,160,320,640,1280,2560];
DOFS = [84,246,810,2898,10914,42306,166530,660738];
DOFS_LST = [246,810,2898,10914,42306,166530,660738,2632194];

% K assembly
figure;
loglog(DOFS(1:5),truss_K_MATLAB,'*-',DOFS,truss_K,'o-',DOFS(1:5),CST_K_MATLAB,'k*-',...
    DOFS,CST_K,'o-',DOFS_LST(1:4),LST_K_MATLAB,'*-',DOFS_LST(1:end-1),LST_K(1:end-1),'o-');
grid on;
xlabel('DOFs'); ylabel('Time [s]');
title('K Assembly Performance');
legend('MATLAB Truss','C++ Truss','MATLAB CST','C++ CST','MATLAB LST','C++ LST');


% K inversion
figure;

loglog(DOFS(1:5),truss_inv_MATLAB,'*-',DOFS,truss_inv,'o-',DOFS(1:5),CST_inv_MATLAB,...
    'k*-',DOFS,CST_inv,'o-',DOFS_LST(1:4),LST_inv_MATLAB,'*-',DOFS_LST(1:end-1),LST_inv(1:end-1),'o-');
grid on;
xlabel('DOFs'); ylabel('Time [s]');
title('K Inversion');
legend('MATLAB Truss','C++ Truss','MATLAB CST','C++ CST','MATLAB LST','C++ LST');
% p = polyfit(log(DOFS(3:5)),log(LST_inv_MATLAB(2:end)),1);
% z = polyval(p,log(DOFS(3:5)));
% hold on;loglog(DOFS(3:5),exp(z))

% Optimisation
% figure;
% loglog(DOFS,truss_optimisation,'o-',DOFS,CST_optimisation,'o-',DOFS_LST,LST_optimisation,'o-');
% grid on;
% xlabel('DOFs'); ylabel('Times [s]');
% title('Optimisation');
% legend('Truss','CST','LST');

figure;
semilogx(NX,truss_optimisation,'o-',NX,CST_optimisation,'o-',NX,LST_optimisation,'ko-');
grid on;
xlabel('n_x'); ylabel('Time [s]');
title('Optimisation');
legend('Truss','CST','LST');