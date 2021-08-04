% Using Lecture_4_2Dtruss by Dr. Sadowski as a template
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimensions of domain and elements in each direction
nx = 80*12;
ny = 4*12;
n_el = nx*2*ny;
Lx = 2000;
Ly = 100;
dx = Lx/nx; % Distance between nodes in x
dy = Ly/ny; % Distance between nodes in y
Nx = nx+1; % Nodes in x direction
Ny = ny+1; % Nodes in y direction
N = Nx*Ny; % Total number of nodes

% Parameters
E = 2e5; % N/m2 - modulus of elasticity of each bar element
nu = 0.3; % Poisson coefficient
t = 10; % Element thickness in mm
loadw = -1e5; % N/m2 - Imposed UDL
% q = loadw.*ones(nx,1);

q_indices = randi([0 1],nx,1);
q_values = loadw.*rand(nx,1);
q = q_values.*q_indices;
weight = -80e3; % N/m3 - Unit weight of steel

% Specifying nodal x-y coordinates
NODES.coords = zeros(N,2);
for i=1:Nx
    for j=1:Ny
        NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
    end
end

NODES.dofs = zeros(N,2);
for i=1:N  
    NODES.dofs(i,:) = [2*i-1, 2*i];
end
% Note that node numbers are identified by their row number
     
% Specifying element nodal connectivity (order does not matter)
ELEMENTS = zeros(n_el,3);
for i=1:ny
    for j=1:nx
        ELEMENTS((i-1)*2*nx+2*j-1,:) = [i+(j-1)*(ny+1), i+ny+1+(j-1)*(ny+1), i+1+(j-1)*(ny+1)];
        ELEMENTS((i-1)*2*nx+2*j,:) = [i+(j-1)*(ny+1)+1, i+ny+1+(j-1)*(ny+1), i+ny+2+(j-1)*(ny+1)];
    end
end

% Degrees of freedom & other parameters
dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*Ny-2,2*Ny:2*N]; % unknown nodal x-y dofs
dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*Ny-1]; % known nodal x-y dofs due to BC (at nodes 1,9)
nodes = size(NODES.coords,1); % no. of nodes i.e. 6
elements = size(ELEMENTS,1); % no. of elements i.e.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the global stiffness matrix
tic;
counter = 1;
F = zeros(2*nodes,1);
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); % identify element node numbers
    n3 = ELEMENTS(EL,3);
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y coordinates
    
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2); % element node 1 - dofs
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2); % element node 2 - dofs
    dof31 = NODES.dofs(n3,1); dof32 = NODES.dofs(n3,2); % element node 3 - dofs
    
    % Nodal force vector
    x21 = x2 - x1; x31 = x3 - x1; % Triangle sides
    y21 = y2 - y1; y31 = y3 - y1;
    A = abs(x21*y31 - x31*y21)/2; % Area of element
    
    F(dof12) = F(dof12) + A*weight*t*1e-9/3; % Weight contribution
    F(dof22) = F(dof22) + A*weight*t*1e-9/3;
    F(dof32) = F(dof32) + A*weight*t*1e-9/3;
    
    if y1 == Ly  && y3 == Ly % UDL contribution
        F(dof12) = F(dof12) + q(counter)*dx*1e-3/2;
        F(dof32) = F(dof32) + q(counter)*dx*1e-3/2;
        counter = counter + 1;
    end
end
disp('Hola ');

coords = reshape(NODES.coords',1,[]);
elem = uint64(reshape(ELEMENTS',1,[]));


% mex CST_K.cpp
tic;
K_sparse = CST_K(E,nu,t,coords,elem);

t1 = toc;

disp(['C++ exec time = ',num2str(t1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SOLVER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specification of submatrices directly, note there is no need to form
% a rearranged K or F explicitly
KRR = K_sparse(dofs_restrained,dofs_restrained);
KRF = K_sparse(dofs_restrained,dofs_free);
KFR = K_sparse(dofs_free,dofs_restrained);
KFF = K_sparse(dofs_free,dofs_free);
fF = F(dofs_free);
% You cannot yet form fR as you do not know what the base reaction is.
% Also, you can even avoid formulating KRR etc. separately and just access
% the requires sub-matrices of K and f directly in what follows.

% Solution for the unknown nodal dofs
uR = zeros(length(dofs_restrained),1); % BC - zero displacement on nodes 1,2,3,4,5,6,7,8,11,12
uF = KFF\(fF - KFR*uR); % 1st matrix equation
U = zeros(2*nodes,1);
U(dofs_restrained) = uR;
U(dofs_free) = uF; % full nodal dof vector

% Solution for the unknown reactions
fR = KRF*uF + KRR*uR; % 2nd matrix equation
F(dofs_restrained) = fR; % full nodal force vector

% Get nodal force vector
% forces = K_sparse*U;

x0 = zeros(nx,1);

tic;
[vars,fval,rel_tol,it] = GN_CST(x0,t,weight,Ly,coords,elem,NODES.dofs,dofs_free,F,1e-4,1e6,0);
t1 = toc;

disp(['Imposed UDL is = ',num2str(vars(1)/1e3),'N/m']);
disp(['Optimisation time = ',num2str(t1)]);