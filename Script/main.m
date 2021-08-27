% Iñigo Basterretxea Jacob
% 01246662
% This file reads the user input parametres from a text file and does
% various things depending on the optins chosen:

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read file
fileID = fopen('input.txt','r');
formatSpec = '%f';
params = fscanf(fileID,formatSpec);

% In case input is incomplete
if length(params) == 20
    % Dimensions of domain and elements in each direction
    Lx = params(1); % Height of domain in mm
    Ly = params(2); % Width of domain in mm
    nx = params(3); % no. of divisions along x-direction
    ny = params(4); % no. of divisions along y-direction
    
    % Parameters
    E = params(5); % N/mm2 - modulus of elasticity of each bar element
    nu = params(7); % Poisson ratio
    loadw = params(8); % N/m - Imposed UDL
    weight = params(9); % N/m3 - Unit weight of material
    load_type = params(10); % Uniform or random set of UDLs
    BC = params(11); % Boundary conditions (1 cantilever or 0 simply-supported)
    element_type = params(12); % 0 truss, 1 CST, 2LST
    stiffness = params(13); % 1 for assembling the stiffness matrix K
    solve_U = params(14); % 1 for solving for the nodal displacements U
    reconstruct_F = params(15); % 1 for load reconstruction
    tol = params(16); % input L2-norm tolerance
    maxIt = params(17); % maximum no. of iterations
    lambda = params(18); % Levenberg–Marquardt Modification coefficient
    plot_U = params(19); % 1 for plotting the elements with(out) nodal displacements
    amp = params(end); % amplification factor for plotting purposes only
    
    if BC == 0 && mod(ny,2) ~= 0 && element_type ~= 2
        msg = 'ABORTED: pure simply-supported conditions (at Neutral Axis) are incompatible with an odd no. of y-divisions (ny) for truss and CST elements';
        error(msg);
    end
    
else
    check = 1;
    while check
        fprintf('1 OR MORE INPUTS MISSING\n');
        default = input('Would you like to go ahead with default parametres (1 = YES)?\n');
        if default
            Lx = 2000;
            Ly = 100;
            nx = 40;
            ny = 2;
            E = 2e5;
            nu = 0.3;
            loadw = -1e4;
            weight = -80e3;
            load_type = 0;
            BC = 1;
            element_type = 1;
            stiffness = 1;
            solve_U = 1;
            reconstruct_F = 1;
            tol = 1e-6;
            maxIt = 1000;
            lambda = 0;
            plot_U = 1;
            amp = 1;
            check = 0;
        else
            msg = 'ABORTED due to lack of input.';
            error(msg);
        end
    end
end

if element_type == 1 || element_type == 2
    t = params(6); % Element thickness in mm
else
    t = params(6)^2*pi; % mm2 - cross-sectional area of each bar element
end

if load_type == 1
    q_indices = randi([0 1],nx,1);
    q_values = loadw.*rand(nx,1);
    q = q_values.*q_indices; % Random vector of upper boundary UDLs
else
    q = loadw.*ones(nx,1); % Uniform vector of upper boundary UDLs
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% MESH %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get mesh (i.e. node, element nodal connectivity, etc.)
[NODES,ELEMENTS,dofs_free,dofs_restrained] = mesh_FE(Lx,Ly,nx,ny,element_type,BC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ASSEMBLE K %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stiffness
    K = K_FE(ELEMENTS,NODES.coords,E,t,nu,element_type);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ASSEMBLE F %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = F_FE(ELEMENTS,NODES.coords,NODES.dofs,Ly,t,q,weight,element_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SOLVER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We obviously cannot solve for U unless we have the stiffness matrix
if solve_U && stiffness 
    [U,F] = solve_FE(K,F,dofs_free,dofs_restrained);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% OPTIMISATION MODULE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if reconstruct_F
    % Package parametres
    options = [tol,maxIt,lambda];
    constants = [Ly,t,weight];
    
    % Initial guess
    x0 = zeros(nx,1);
    [x_sol,fval,it] = load_reconstruction(x0,F,ELEMENTS,NODES.coords,NODES.dofs,dofs_free,options,constants,element_type);

    % Calculate maximum relative error
    error = abs(x_sol - q)./abs(q).*100;
    error(isinf(error)) = 0;
    max_error = max(error);
    if max_error < 1e-6
        disp(['CORRECT: maximum relative error is = ',num2str(max_error),'%']);
    else
        disp(['ERROR: maximum relative error is larger than = ',num2str(1e-6),'%']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% POST-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_U
    p = plot_FE(ELEMENTS,NODES.coords,NODES.dofs,U,amp,element_type,solve_U);
end