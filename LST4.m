% Using Lecture_4_2Dtruss by Dr. Sadowski as a template
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
E = 2e5; % N/m2 - modulus of elasticity of each bar element
nu = 0.3; % Poisson coefficient
t = 10; % Element thickness in mm

% Dimensions of domain and elements in each direction
nx = 20;
ny = 1;
n_el = nx*2*ny;
Lx = 2000;
Ly = 100;
dx = Lx/(2*nx); % Distance between nodes in x
dy = Ly/(2*ny); % Distance between nodes in y
Nx = 2*nx+1; % Nodes in x direction
Ny = 2*ny+1; % Nodes in y direction
N = Nx*Ny; % Total number of nodes

weight = -80e3; % N/m3 - Unit weight of steel
loadw = -1e5; % N/m2 - Imposed UDL
% q = loadw.*ones(nx,1);

q_indices = randi([0 1],nx,1);
q_values = loadw.*rand(nx,1);
q = q_values.*q_indices;

% Specifying nodal x-y coordinates
NODES.coords = zeros(N,2);
for i=1:Ny
    for j=1:Nx
        NODES.coords((i-1)*Nx+j,:) = [(j-1)*dx,dy*(i-1)];
    end
end

NODES.dofs = zeros(N,2);
for i=1:N
    
    NODES.dofs(i,:) = [2*i-1, 2*i];
end
% Note that node numbers are identified by their row number
     
% Specifying element nodal connectivity (order does not matter)
ELEMENTS = zeros(n_el,6);
for i=1:n_el
    if mod(floor(i/nx),2) == 0 && mod(i,nx)~= 0
        ELEMENTS(i,:) = [2*mod(i,nx)-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0)+1, 2*Nx+2*mod(i,nx)-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), Nx+2*mod(i,nx)+Nx*floor(i/nx), Nx+2*mod(i,nx)-1+Nx*floor(i/nx)];
    elseif mod(i,nx)== 0 && mod(floor(i/nx),2) ~= 0
        ELEMENTS(i,:) = [Nx-2+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), Nx+2*Nx-2+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), Nx-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), -1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), -2+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0)]; %Nx+2*Nx-2+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0)
    elseif mod(i,nx)== 0 && mod(floor(i/nx),2) == 0
        ELEMENTS(i,:) = [2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)+Nx, Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-2, Nx+2*Nx*(max(floor(i/nx)-floor(i/(2*nx)),1)-1), Nx-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1), Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx-1, Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx];
    else
        ELEMENTS(i,:) = [2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)+2*mod(i,nx)+1, 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-1, 2*mod(i,nx)+1+2*Nx*(max(floor(i/nx)-floor(i/(2*nx)),1)-1), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx, 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx+1];
    end  
end

% Degrees of freedom & other parameters
dofs_free = [1:N-Nx,N-Nx+3:N+Nx-1,N+Nx+1:2*N]; % unknown nodal x-y dofs
dofs_restrained = [N-Nx+1,N-Nx+2,N+Nx]; % known nodal x-y dofs due to BC (at nodes 1,9)
nodes = size(NODES.coords,1); % no. of nodes i.e. 52
elements = size(ELEMENTS,1); % no. of elements i.e. 16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the global nodal force vector
F = zeros(2*nodes,1); % initialising en empty a 2Nx1 column vector for convenience
% % Constructing the global stiffness matrix
% K = zeros(2*nodes); % initialising en empty 66x66 matrix; 33 nodes at 2 dofs/node
count = 1;
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); % identify element node numbers
    n4 = ELEMENTS(EL,4); n5 = ELEMENTS(EL,5); n6 = ELEMENTS(EL,6); % identify element node numbers
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y coordinates
    x4 = NODES.coords(n4,1); y4 = NODES.coords(n4,2); % element node 1 - x,y coordinates
    x5 = NODES.coords(n5,1); y5 = NODES.coords(n5,2); % element node 2 - x,y coordinates
    x6 = NODES.coords(n6,1); y6 = NODES.coords(n6,2); % element node 3 - x,y coordinates
    
%     b = abs(y2-y1); h = abs(x3-x1); % base and height of triangular element
%     A = 0.5*abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)); % area of triangle
    
    
%     dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2);  % element node 1 - dofs
%     dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2);  % element node 2 - dofs
%     dof31 = NODES.dofs(n3,1); dof32 = NODES.dofs(n3,2);  % element node 3 - dofs
%     dof41 = NODES.dofs(n4,1); dof42 = NODES.dofs(n4,2);  % element node 4 - dofs
%     dof51 = NODES.dofs(n5,1); dof52 = NODES.dofs(n5,2);  % element node 5 - dofs
%     dof61 = NODES.dofs(n6,1); dof62 = NODES.dofs(n6,2);  % element node 6 - dofs
%         
%     constants = E*t/(1-nu.^2); % Constants over element
%     
%     ke = k_mat(x1,x2,x3,y1,y2,y3,nu);
% 
%     ke = constants*ke; % element stiffness
%         
%     
%     % Updating global stiffness matrix [K] coefficients 
%     % Note that for each element you still have to 'think locally'
%     % Row 1 - element dof11
%     K(dof11,dof11) = K(dof11,dof11) + ke(1,1); % Col 1 - element dof11
%     K(dof11,dof12) = K(dof11,dof12) + ke(1,2); % Col 2 - element dof12
%     K(dof11,dof21) = K(dof11,dof21) + ke(1,3); % Col 3 - element dof21
%     K(dof11,dof22) = K(dof11,dof22) + ke(1,4); % Col 4 - element dof22
%     K(dof11,dof31) = K(dof11,dof31) + ke(1,5); % Col 5 - element dof31
%     K(dof11,dof32) = K(dof11,dof32) + ke(1,6); % Col 6 - element dof32
%     K(dof11,dof41) = K(dof11,dof41) + ke(1,7); % Col 7 - element dof41
%     K(dof11,dof42) = K(dof11,dof42) + ke(1,8); % Col 8 - element dof42
%     K(dof11,dof51) = K(dof11,dof51) + ke(1,9); % Col 9 - element dof51
%     K(dof11,dof52) = K(dof11,dof52) + ke(1,10); % Col 10 - element dof52
%     K(dof11,dof61) = K(dof11,dof61) + ke(1,11); % Col 11 - element dof61
%     K(dof11,dof62) = K(dof11,dof62) + ke(1,12); % Col 12 - element dof62
%     
%     % Row 2 - element dof12
%     K(dof12,dof11) = K(dof12,dof11) + ke(2,1); % Col 1 - element dof11
%     K(dof12,dof12) = K(dof12,dof12) + ke(2,2); % Col 2 - element dof12
%     K(dof12,dof21) = K(dof12,dof21) + ke(2,3); % Col 3 - element dof21
%     K(dof12,dof22) = K(dof12,dof22) + ke(2,4); % Col 4 - element dof22
%     K(dof12,dof31) = K(dof12,dof31) + ke(2,5); % Col 5 - element dof31
%     K(dof12,dof32) = K(dof12,dof32) + ke(2,6); % Col 6 - element dof32
%     K(dof12,dof41) = K(dof12,dof41) + ke(2,7); % Col 7 - element dof41
%     K(dof12,dof42) = K(dof12,dof42) + ke(2,8); % Col 8 - element dof42
%     K(dof12,dof51) = K(dof12,dof51) + ke(2,9); % Col 9 - element dof51
%     K(dof12,dof52) = K(dof12,dof52) + ke(2,10); % Col 10 - element dof52
%     K(dof12,dof61) = K(dof12,dof61) + ke(2,11); % Col 11 - element dof61
%     K(dof12,dof62) = K(dof12,dof62) + ke(2,12); % Col 12 - element dof62
%     
%     % Row 3 - element dof21
%     K(dof21,dof11) = K(dof21,dof11) + ke(3,1); % Col 1 - element dof11
%     K(dof21,dof12) = K(dof21,dof12) + ke(3,2); % Col 2 - element dof12
%     K(dof21,dof21) = K(dof21,dof21) + ke(3,3); % Col 3 - element dof21
%     K(dof21,dof22) = K(dof21,dof22) + ke(3,4); % Col 4 - element dof22
%     K(dof21,dof31) = K(dof21,dof31) + ke(3,5); % Col 5 - element dof31
%     K(dof21,dof32) = K(dof21,dof32) + ke(3,6); % Col 6 - element dof32
%     K(dof21,dof41) = K(dof21,dof41) + ke(3,7); % Col 7 - element dof41
%     K(dof21,dof42) = K(dof21,dof42) + ke(3,8); % Col 8 - element dof42
%     K(dof21,dof51) = K(dof21,dof51) + ke(3,9); % Col 9 - element dof51
%     K(dof21,dof52) = K(dof21,dof52) + ke(3,10); % Col 10 - element dof52
%     K(dof21,dof61) = K(dof21,dof61) + ke(3,11); % Col 11 - element dof61
%     K(dof21,dof62) = K(dof21,dof62) + ke(3,12); % Col 12 - element dof62
%     
%     % Row 4 - element dof22
%     K(dof22,dof11) = K(dof22,dof11) + ke(4,1); % Col 1 - element dof11
%     K(dof22,dof12) = K(dof22,dof12) + ke(4,2); % Col 2 - element dof12
%     K(dof22,dof21) = K(dof22,dof21) + ke(4,3); % Col 3 - element dof21
%     K(dof22,dof22) = K(dof22,dof22) + ke(4,4); % Col 4 - element dof22
%     K(dof22,dof31) = K(dof22,dof31) + ke(4,5); % Col 5 - element dof31
%     K(dof22,dof32) = K(dof22,dof32) + ke(4,6); % Col 6 - element dof32
%     K(dof22,dof41) = K(dof22,dof41) + ke(4,7); % Col 7 - element dof41
%     K(dof22,dof42) = K(dof22,dof42) + ke(4,8); % Col 8 - element dof42
%     K(dof22,dof51) = K(dof22,dof51) + ke(4,9); % Col 9 - element dof51
%     K(dof22,dof52) = K(dof22,dof52) + ke(4,10); % Col 10 - element dof52
%     K(dof22,dof61) = K(dof22,dof61) + ke(4,11); % Col 11 - element dof61
%     K(dof22,dof62) = K(dof22,dof62) + ke(4,12); % Col 12 - element dof62
%     
%     % Row 5 - element dof31
%     K(dof31,dof11) = K(dof31,dof11) + ke(5,1); % Col 1 - element dof11
%     K(dof31,dof12) = K(dof31,dof12) + ke(5,2); % Col 2 - element dof12
%     K(dof31,dof21) = K(dof31,dof21) + ke(5,3); % Col 3 - element dof21
%     K(dof31,dof22) = K(dof31,dof22) + ke(5,4); % Col 4 - element dof22
%     K(dof31,dof31) = K(dof31,dof31) + ke(5,5); % Col 5 - element dof31
%     K(dof31,dof32) = K(dof31,dof32) + ke(5,6); % Col 6 - element dof32
%     K(dof31,dof41) = K(dof31,dof41) + ke(5,7); % Col 7 - element dof41
%     K(dof31,dof42) = K(dof31,dof42) + ke(5,8); % Col 8 - element dof42
%     K(dof31,dof51) = K(dof31,dof51) + ke(5,9); % Col 9 - element dof51
%     K(dof31,dof52) = K(dof31,dof52) + ke(5,10); % Col 10 - element dof52
%     K(dof31,dof61) = K(dof31,dof61) + ke(5,11); % Col 11 - element dof61
%     K(dof31,dof62) = K(dof31,dof62) + ke(5,12); % Col 12 - element dof62
%     
%     % Row 6 - element dof32
%     K(dof32,dof11) = K(dof32,dof11) + ke(6,1); % Col 1 - element dof11
%     K(dof32,dof12) = K(dof32,dof12) + ke(6,2); % Col 2 - element dof12
%     K(dof32,dof21) = K(dof32,dof21) + ke(6,3); % Col 3 - element dof21
%     K(dof32,dof22) = K(dof32,dof22) + ke(6,4); % Col 4 - element dof22
%     K(dof32,dof31) = K(dof32,dof31) + ke(6,5); % Col 5 - element dof31
%     K(dof32,dof32) = K(dof32,dof32) + ke(6,6); % Col 6 - element dof32
%     K(dof32,dof41) = K(dof32,dof41) + ke(6,7); % Col 7 - element dof41
%     K(dof32,dof42) = K(dof32,dof42) + ke(6,8); % Col 8 - element dof42
%     K(dof32,dof51) = K(dof32,dof51) + ke(6,9); % Col 9 - element dof51
%     K(dof32,dof52) = K(dof32,dof52) + ke(6,10); % Col 10 - element dof52
%     K(dof32,dof61) = K(dof32,dof61) + ke(6,11); % Col 11 - element dof61
%     K(dof32,dof62) = K(dof32,dof62) + ke(6,12); % Col 12 - element dof62
%     
%     % Row 7 - element dof41
%     K(dof41,dof11) = K(dof41,dof11) + ke(7,1); % Col 1 - element dof11
%     K(dof41,dof12) = K(dof41,dof12) + ke(7,2); % Col 2 - element dof12
%     K(dof41,dof21) = K(dof41,dof21) + ke(7,3); % Col 3 - element dof21
%     K(dof41,dof22) = K(dof41,dof22) + ke(7,4); % Col 4 - element dof22
%     K(dof41,dof31) = K(dof41,dof31) + ke(7,5); % Col 5 - element dof31
%     K(dof41,dof32) = K(dof41,dof32) + ke(7,6); % Col 6 - element dof32
%     K(dof41,dof41) = K(dof41,dof41) + ke(7,7); % Col 7 - element dof41
%     K(dof41,dof42) = K(dof41,dof42) + ke(7,8); % Col 8 - element dof42
%     K(dof41,dof51) = K(dof41,dof51) + ke(7,9); % Col 9 - element dof51
%     K(dof41,dof52) = K(dof41,dof52) + ke(7,10); % Col 10 - element dof52
%     K(dof41,dof61) = K(dof41,dof61) + ke(7,11); % Col 11 - element dof61
%     K(dof41,dof62) = K(dof41,dof62) + ke(7,12); % Col 12 - element dof62
%     
%     % Row 8 - element dof42
%     K(dof42,dof11) = K(dof42,dof11) + ke(8,1); % Col 1 - element dof11
%     K(dof42,dof12) = K(dof42,dof12) + ke(8,2); % Col 2 - element dof12
%     K(dof42,dof21) = K(dof42,dof21) + ke(8,3); % Col 3 - element dof21
%     K(dof42,dof22) = K(dof42,dof22) + ke(8,4); % Col 4 - element dof22
%     K(dof42,dof31) = K(dof42,dof31) + ke(8,5); % Col 5 - element dof31
%     K(dof42,dof32) = K(dof42,dof32) + ke(8,6); % Col 6 - element dof32
%     K(dof42,dof41) = K(dof42,dof41) + ke(8,7); % Col 7 - element dof41
%     K(dof42,dof42) = K(dof42,dof42) + ke(8,8); % Col 8 - element dof42
%     K(dof42,dof51) = K(dof42,dof51) + ke(8,9); % Col 9 - element dof51
%     K(dof42,dof52) = K(dof42,dof52) + ke(8,10); % Col 10 - element dof52
%     K(dof42,dof61) = K(dof42,dof61) + ke(8,11); % Col 11 - element dof61
%     K(dof42,dof62) = K(dof42,dof62) + ke(8,12); % Col 12 - element dof62
%     
%     % Row 9 - element dof51
%     K(dof51,dof11) = K(dof51,dof11) + ke(9,1); % Col 1 - element dof11
%     K(dof51,dof12) = K(dof51,dof12) + ke(9,2); % Col 2 - element dof12
%     K(dof51,dof21) = K(dof51,dof21) + ke(9,3); % Col 3 - element dof21
%     K(dof51,dof22) = K(dof51,dof22) + ke(9,4); % Col 4 - element dof22
%     K(dof51,dof31) = K(dof51,dof31) + ke(9,5); % Col 5 - element dof31
%     K(dof51,dof32) = K(dof51,dof32) + ke(9,6); % Col 6 - element dof32
%     K(dof51,dof41) = K(dof51,dof41) + ke(9,7); % Col 7 - element dof41
%     K(dof51,dof42) = K(dof51,dof42) + ke(9,8); % Col 8 - element dof42
%     K(dof51,dof51) = K(dof51,dof51) + ke(9,9); % Col 9 - element dof51
%     K(dof51,dof52) = K(dof51,dof52) + ke(9,10); % Col 10 - element dof52
%     K(dof51,dof61) = K(dof51,dof61) + ke(9,11); % Col 11 - element dof61
%     K(dof51,dof62) = K(dof51,dof62) + ke(9,12); % Col 12 - element dof62
%     
%     % Row 10 - element dof52
%     K(dof52,dof11) = K(dof52,dof11) + ke(10,1); % Col 1 - element dof11
%     K(dof52,dof12) = K(dof52,dof12) + ke(10,2); % Col 2 - element dof12
%     K(dof52,dof21) = K(dof52,dof21) + ke(10,3); % Col 3 - element dof21
%     K(dof52,dof22) = K(dof52,dof22) + ke(10,4); % Col 4 - element dof22
%     K(dof52,dof31) = K(dof52,dof31) + ke(10,5); % Col 5 - element dof31
%     K(dof52,dof32) = K(dof52,dof32) + ke(10,6); % Col 6 - element dof32
%     K(dof52,dof41) = K(dof52,dof41) + ke(10,7); % Col 7 - element dof41
%     K(dof52,dof42) = K(dof52,dof42) + ke(10,8); % Col 8 - element dof42
%     K(dof52,dof51) = K(dof52,dof51) + ke(10,9); % Col 9 - element dof51
%     K(dof52,dof52) = K(dof52,dof52) + ke(10,10); % Col 10 - element dof52
%     K(dof52,dof61) = K(dof52,dof61) + ke(10,11); % Col 11 - element dof61
%     K(dof52,dof62) = K(dof52,dof62) + ke(10,12); % Col 12 - element dof62
%     
%     % Row 11 - element dof61
%     K(dof61,dof11) = K(dof61,dof11) + ke(11,1); % Col 1 - element dof11
%     K(dof61,dof12) = K(dof61,dof12) + ke(11,2); % Col 2 - element dof12
%     K(dof61,dof21) = K(dof61,dof21) + ke(11,3); % Col 3 - element dof21
%     K(dof61,dof22) = K(dof61,dof22) + ke(11,4); % Col 4 - element dof22
%     K(dof61,dof31) = K(dof61,dof31) + ke(11,5); % Col 5 - element dof31
%     K(dof61,dof32) = K(dof61,dof32) + ke(11,6); % Col 6 - element dof32
%     K(dof61,dof41) = K(dof61,dof41) + ke(11,7); % Col 7 - element dof41
%     K(dof61,dof42) = K(dof61,dof42) + ke(11,8); % Col 8 - element dof42
%     K(dof61,dof51) = K(dof61,dof51) + ke(11,9); % Col 9 - element dof51
%     K(dof61,dof52) = K(dof61,dof52) + ke(11,10); % Col 10 - element dof52
%     K(dof61,dof61) = K(dof61,dof61) + ke(11,11); % Col 11 - element dof61
%     K(dof61,dof62) = K(dof61,dof62) + ke(11,12); % Col 12 - element dof62
%     
%     % Row 12 - element dof62
%     K(dof62,dof11) = K(dof62,dof11) + ke(12,1); % Col 1 - element dof11
%     K(dof62,dof12) = K(dof62,dof12) + ke(12,2); % Col 2 - element dof12
%     K(dof62,dof21) = K(dof62,dof21) + ke(12,3); % Col 3 - element dof21
%     K(dof62,dof22) = K(dof62,dof22) + ke(12,4); % Col 4 - element dof22
%     K(dof62,dof31) = K(dof62,dof31) + ke(12,5); % Col 5 - element dof31
%     K(dof62,dof32) = K(dof62,dof32) + ke(12,6); % Col 6 - element dof32
%     K(dof62,dof41) = K(dof62,dof41) + ke(12,7); % Col 7 - element dof41
%     K(dof62,dof42) = K(dof62,dof42) + ke(12,8); % Col 8 - element dof42
%     K(dof62,dof51) = K(dof62,dof51) + ke(12,9); % Col 9 - element dof51
%     K(dof62,dof52) = K(dof62,dof52) + ke(12,10); % Col 10 - element dof52
%     K(dof62,dof61) = K(dof62,dof61) + ke(12,11); % Col 11 - element dof61
%     K(dof62,dof62) = K(dof62,dof62) + ke(12,12); % Col 12 - element dof62
    
    dof12 = NODES.dofs(n1,2); % element node 1 - dofs
    dof22 = NODES.dofs(n2,2); % element node 2 - dofs
    dof42 = NODES.dofs(n4,2); % element node 3 - dofs
    dof52 = NODES.dofs(n5,2); % element node 3 - dofs
    dof62 = NODES.dofs(n6,2); % element node 3 - dofs
    
    x21 = x2 - x1; x31 = x3 - x1; % Triangle sides
    y21 = y2 - y1; y31 = y3 - y1;
    A = abs(x21*y31 - x31*y21)/2; % Area of element
    
    x41 = x4 - x1;
    
    % Weight contribution. No vertex loads!
    F(dof42) = F(dof42) + A*weight*t*1e-9/3; 
    F(dof52) = F(dof52) + A*weight*t*1e-9/3;
    F(dof62) = F(dof62) + A*weight*t*1e-9/3;
    
    if y1 == Ly  && y2 == Ly && y4 == Ly % UDL contribution
        F(dof12) = F(dof12) + q(count)*abs(x41)*1e-3/2;
        F(dof22) = F(dof22) + q(count)*abs(x41)*1e-3/2;
        % Double contribution for midpoint along edge
        F(dof42) = F(dof42) + q(count)*abs(x41)*1e-3; 
        count = count + 1;
        
    end
    
    
end

coords = reshape(NODES.coords',1,[]);
elem = uint32(reshape(ELEMENTS',1,[]));

tic;

[k_values,k_rows,k_cols]  = LST_K(E,nu,t,coords,elem);
k_rows = double(k_rows) + ones(size(k_rows));
k_cols = double(k_cols) + ones(size(k_cols));
K = sparse(k_rows,k_cols,k_values,2*nodes,2*nodes);

t1 = toc;

disp(['C++ exec time = ',num2str(t1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SOLVER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specification of submatrices directly, note there is no need to form
% a rearranged K or F explicitly
KRR = K(dofs_restrained,dofs_restrained);
KRF = K(dofs_restrained,dofs_free);
KFR = K(dofs_free,dofs_restrained);
KFF = K(dofs_free,dofs_free);
fF = F(dofs_free);
% You cannot yet form fR as you do not know what the base reaction is.
% Also, you can even avoid formulating KRR etc. separately and just access
% the requires sub-matrices of K and f directly in what follows.

% Solution for the unknown nodal dofs
uR = zeros(length(dofs_restrained),1); % BC - zero displacement on nodes 1,5
uF = KFF\(fF - KFR*uR); % 1st matrix equation
U = zeros(2*nodes,1);
U(dofs_restrained) = uR;
U(dofs_free) = uF; % full nodal dof vector

% Solution for the unknown reactions
fR = KRF*uF + KRR*uR; % 2nd matrix equation
F(dofs_restrained) = fR; % full nodal force vector


x0 = zeros(nx,1);

tic;
[vars,fval,rel_tol,it] = GN_LST(x0,t,weight,Ly,coords,elem,NODES.dofs,dofs_free,F,1e-4,1e6,0);
t1 = toc;

disp(['Imposed UDL is = ',num2str(vars(1)/1e3),'N/m']);
disp(['Optimisation time = ',num2str(t1)]);