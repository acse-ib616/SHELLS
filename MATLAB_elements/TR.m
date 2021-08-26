% Iñigo Basterretxea Jacob 01246662

% Script that implements the Finite Element Method for a rectangular domain
% with TR elements
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
E = 2e5; % N/mm2 - modulus of elasticity of each bar element
nu = 0.3; % Poisson ratio
t = 10; % Element thickness in mm
weight = -80e3; % N/m3 - Unit weight of steel
loadw = -1e4; % N/m2 - Imposed UDL

% Dimensions of domain and elements in each direction
nx = 40; % no. of divisions along x-direction
ny = 20; % no. of divisions along y-direction
Lx = 2000; % Height of domain in mm
Ly = 100; % Width of domain in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx_LST = Lx/(2*nx); % Distance between nodes in x (LST)
dy_LST = Ly/(2*ny); % Distance between nodes in y (LST)
Nx_LST = 2*nx+1; % no. of nodes in x direction (LST)
Ny_LST = 2*ny+1; % no. of nodes in y direction (LST)
N_LST = Nx_LST*Ny_LST; % no. of nodes number of nodes (LST)

elements = nx*2*ny; % no. of elements
Nx = nx+1; % Nodes in x direction (TR)
Ny = ny+1; % Nodes in y direction (TR)
nodes = Nx*Ny; % no. of nodes (TR)
dx = Lx/nx; % Distance between nodes in x
dy = Ly/ny; % Distance between nodes in y

% Note that node numbers are identified by their row number
% Specifying nodal x-y coordinates. 1st column is x & 2nd is y
NODES.coords = zeros(nodes,2);
for i=1:Nx
    for j=1:Ny
        NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
    end
end

% Specifying  DOFs in 3s per node. 1st column are horizontal DOFs, 2nd
% column are vertical DOFs & 3rd column are rotational DOFs
NODES.dofs = zeros(nodes,3);
for i=1:nodes
    NODES.dofs(i,:) = [3*i-2, 3*i-1, 3*i];
end

% Specifying element nodal connectivity in rows of 3 nodes per element
% (order of the elements does not matter)
ELEMENTS = zeros(elements,3);
for i=1:ny
    for j=1:nx
        ELEMENTS((i-1)*2*nx+2*j-1,:) = [i+(j-1)*(ny+1), i+ny+1+(j-1)*(ny+1), i+1+(j-1)*(ny+1)];
        ELEMENTS((i-1)*2*nx+2*j,:) = [i+(j-1)*(ny+1)+1, i+ny+1+(j-1)*(ny+1), i+ny+2+(j-1)*(ny+1)];
    end
end


% Degrees of freedom & other parameters
% unknown nodal x-y dofs
dofs_free = [1:(floor(ny/2)+1)*3-3,(floor(ny/2)+1)*3+1:3*Ny-3,3*Ny-1:3*nodes];
% known nodal x-y dofs due to BC
dofs_restrained = [(floor(ny/2)+1)*3-2,(floor(ny/2)+1)*3-1,(floor(ny/2)+1)*3,3*Ny-2]; 

counter = 1; % Initialise q-vector counter
q = loadw.*ones(nx,1); % Vector of upper boundary UDLs
F = zeros(3*nodes,1); % Initialise nodal forces vector

K = zeros(3*nodes); % Initialise global stiffness matrix

% Constructing the global stiffness matrix and the nodal forces vector
for EL = 1:elements % loop through all elements & build stiffness matrix
    
    % identify element node numbers
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); 
    
    % identify node coordinates
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2);
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2);
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2);
    
    % identify element DOFs
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2); dof13 = NODES.dofs(n1,3);
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2); dof23 = NODES.dofs(n2,3);
    dof31 = NODES.dofs(n3,1); dof32 = NODES.dofs(n3,2); dof33 = NODES.dofs(n3,3);
    
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE F %%%
    %%%%%%%%%%%%%%%%%%
    % Triangle sides
    x21 = x2 - x1; x31 = x3 - x1; y21 = y2 - y1; y31 = y3 - y1;
    
    % Area of element
    A = abs(x21*y31 - x31*y21)/2; 
    
    % Weight contribution at all elements
    F(dof12) = F(dof12) + A*weight*t*1e-9/3;
    F(dof22) = F(dof22) + A*weight*t*1e-9/3;
    F(dof32) = F(dof32) + A*weight*t*1e-9/3;
    
    % UDL contribution at top boundary elements
    if y1 == Ly  && y3 == Ly
        F(dof12) = F(dof12) + q(counter)*dx*1e-3/2;
        F(dof32) = F(dof32) + q(counter)*dx*1e-3/2;
        counter = counter + 1;
    end
        
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE K %%%
    %%%%%%%%%%%%%%%%%%
    
    % LST element stiffness matrix
    constants = E*t/(1-nu.^2);
    ke_LST = k_mat_LST(x1,x2,x3,y1,y2,y3,nu);
    ke_LST = constants*ke_LST;
    
    % Transformation from LST to TR
    % Triangle sides
    y32 = y3-y2; x32 = x3-x2; y13 = y1-y3; x13 = x1-x3;
    
    % Transformation matrix
    T = [1 0 0 0 0 0 0 0 0;
         0 1 0 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0;
         0 0 0 0 1 0 0 0 0;
         0 0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 1 0;
         0.5 0 -y21/8 0.5 0 y21/8 0 0 0;
         0 0.5 -x21/8 0 0.5 x21/8 0 0 0;
         0 0 0 0.5 0 -y32/8 0.5 0 y32/8;
         0 0 0 0 0.5 -x32/8 0 0.5 x32/8;
         0.5 0 y13/8 0 0 0 0.5 0 -y13/8;
         0 0.5 x13/8 0 0 0 0 0.5 -x13/8];
     
     ke = T.'*ke_LST*T;
    
    % Updating global stiffness matrix [K] coefficients
    % Row 1 - element dof11
    K(dof11,dof11) = K(dof11,dof11) + ke(1,1); % Col 1 - element dof11
    K(dof11,dof12) = K(dof11,dof12) + ke(1,2); % Col 2 - element dof12
    K(dof11,dof13) = K(dof11,dof13) + ke(1,3); % Col 3 - element dof13
    K(dof11,dof21) = K(dof11,dof21) + ke(1,4); % Col 4 - element dof21
    K(dof11,dof22) = K(dof11,dof22) + ke(1,5); % Col 5 - element dof22
    K(dof11,dof23) = K(dof11,dof23) + ke(1,6); % Col 6 - element dof23
    K(dof11,dof31) = K(dof11,dof31) + ke(1,7); % Col 7 - element dof31
    K(dof11,dof32) = K(dof11,dof32) + ke(1,8); % Col 8 - element dof32
    K(dof11,dof33) = K(dof11,dof33) + ke(1,9); % Col 9 - element dof33
    
    % Row 2 - element dof12
    K(dof12,dof11) = K(dof12,dof11) + ke(2,1); % Col 1 - element dof11
    K(dof12,dof12) = K(dof12,dof12) + ke(2,2); % Col 2 - element dof12
    K(dof12,dof13) = K(dof12,dof13) + ke(2,3); % Col 3 - element dof13
    K(dof12,dof21) = K(dof12,dof21) + ke(2,4); % Col 4 - element dof21
    K(dof12,dof22) = K(dof12,dof22) + ke(2,5); % Col 5 - element dof22
    K(dof12,dof23) = K(dof12,dof23) + ke(2,6); % Col 6 - element dof23
    K(dof12,dof31) = K(dof12,dof31) + ke(2,7); % Col 7 - element dof31
    K(dof12,dof32) = K(dof12,dof32) + ke(2,8); % Col 8 - element dof32
    K(dof12,dof33) = K(dof12,dof33) + ke(2,9); % Col 9 - element dof33
    
    % Row 3 - element dof13
    K(dof13,dof11) = K(dof13,dof11) + ke(3,1); % Col 1 - element dof11
    K(dof13,dof12) = K(dof13,dof12) + ke(3,2); % Col 2 - element dof12
    K(dof13,dof13) = K(dof13,dof13) + ke(3,3); % Col 3 - element dof13
    K(dof13,dof21) = K(dof13,dof21) + ke(3,4); % Col 4 - element dof21
    K(dof13,dof22) = K(dof13,dof22) + ke(3,5); % Col 5 - element dof22
    K(dof13,dof23) = K(dof13,dof23) + ke(3,6); % Col 6 - element dof23
    K(dof13,dof31) = K(dof13,dof31) + ke(3,7); % Col 7 - element dof31
    K(dof13,dof32) = K(dof13,dof32) + ke(3,8); % Col 8 - element dof32
    K(dof13,dof33) = K(dof13,dof33) + ke(3,9); % Col 9 - element dof33
    
    % Row 4 - element dof21
    K(dof21,dof11) = K(dof21,dof11) + ke(4,1); % Col 1 - element dof11
    K(dof21,dof12) = K(dof21,dof12) + ke(4,2); % Col 2 - element dof12
    K(dof21,dof13) = K(dof21,dof13) + ke(4,3); % Col 3 - element dof13
    K(dof21,dof21) = K(dof21,dof21) + ke(4,4); % Col 4 - element dof21
    K(dof21,dof22) = K(dof21,dof22) + ke(4,5); % Col 5 - element dof22
    K(dof21,dof23) = K(dof21,dof23) + ke(4,6); % Col 6 - element dof23
    K(dof21,dof31) = K(dof21,dof31) + ke(4,7); % Col 7 - element dof31
    K(dof21,dof32) = K(dof21,dof32) + ke(4,8); % Col 8 - element dof32
    K(dof21,dof33) = K(dof21,dof33) + ke(4,9); % Col 9 - element dof33
    
    % Row 5 - element dof22
    K(dof22,dof11) = K(dof22,dof11) + ke(5,1); % Col 1 - element dof11
    K(dof22,dof12) = K(dof22,dof12) + ke(5,2); % Col 2 - element dof12
    K(dof22,dof13) = K(dof22,dof13) + ke(5,3); % Col 3 - element dof13
    K(dof22,dof21) = K(dof22,dof21) + ke(5,4); % Col 4 - element dof21
    K(dof22,dof22) = K(dof22,dof22) + ke(5,5); % Col 5 - element dof22
    K(dof22,dof23) = K(dof22,dof23) + ke(5,6); % Col 6 - element dof23
    K(dof22,dof31) = K(dof22,dof31) + ke(5,7); % Col 7 - element dof31
    K(dof22,dof32) = K(dof22,dof32) + ke(5,8); % Col 8 - element dof32
    K(dof22,dof33) = K(dof22,dof33) + ke(5,9); % Col 9 - element dof33
    
    % Row 6 - element dof23
    K(dof23,dof11) = K(dof23,dof11) + ke(6,1); % Col 1 - element dof11
    K(dof23,dof12) = K(dof23,dof12) + ke(6,2); % Col 2 - element dof12
    K(dof23,dof13) = K(dof23,dof13) + ke(6,3); % Col 3 - element dof13
    K(dof23,dof21) = K(dof23,dof21) + ke(6,4); % Col 4 - element dof21
    K(dof23,dof22) = K(dof23,dof22) + ke(6,5); % Col 5 - element dof22
    K(dof23,dof23) = K(dof23,dof23) + ke(6,6); % Col 6 - element dof23
    K(dof23,dof31) = K(dof23,dof31) + ke(6,7); % Col 7 - element dof31
    K(dof23,dof32) = K(dof23,dof32) + ke(6,8); % Col 8 - element dof32
    K(dof23,dof33) = K(dof23,dof33) + ke(6,9); % Col 9 - element dof33
    
    % Row 7 - element dof31
    K(dof31,dof11) = K(dof31,dof11) + ke(7,1); % Col 1 - element dof11
    K(dof31,dof12) = K(dof31,dof12) + ke(7,2); % Col 2 - element dof12
    K(dof31,dof13) = K(dof31,dof13) + ke(7,3); % Col 3 - element dof13
    K(dof31,dof21) = K(dof31,dof21) + ke(7,4); % Col 4 - element dof21
    K(dof31,dof22) = K(dof31,dof22) + ke(7,5); % Col 5 - element dof22
    K(dof31,dof23) = K(dof31,dof23) + ke(7,6); % Col 6 - element dof23
    K(dof31,dof31) = K(dof31,dof31) + ke(7,7); % Col 7 - element dof31
    K(dof31,dof32) = K(dof31,dof32) + ke(7,8); % Col 8 - element dof32
    K(dof31,dof33) = K(dof31,dof33) + ke(7,9); % Col 9 - element dof33
    
    % Row 8 - element dof32
    K(dof32,dof11) = K(dof32,dof11) + ke(8,1); % Col 1 - element dof11
    K(dof32,dof12) = K(dof32,dof12) + ke(8,2); % Col 2 - element dof12
    K(dof32,dof13) = K(dof32,dof13) + ke(8,3); % Col 3 - element dof13
    K(dof32,dof21) = K(dof32,dof21) + ke(8,4); % Col 4 - element dof21
    K(dof32,dof22) = K(dof32,dof22) + ke(8,5); % Col 5 - element dof22
    K(dof32,dof23) = K(dof32,dof23) + ke(8,6); % Col 6 - element dof23
    K(dof32,dof31) = K(dof32,dof31) + ke(8,7); % Col 7 - element dof31
    K(dof32,dof32) = K(dof32,dof32) + ke(8,8); % Col 8 - element dof32
    K(dof32,dof33) = K(dof32,dof33) + ke(8,9); % Col 9 - element dof33
     
    % Row 9 - element dof33
    K(dof33,dof11) = K(dof33,dof11) + ke(9,1); % Col 1 - element dof11
    K(dof33,dof12) = K(dof33,dof12) + ke(9,2); % Col 2 - element dof12
    K(dof33,dof13) = K(dof33,dof13) + ke(9,3); % Col 3 - element dof13
    K(dof33,dof21) = K(dof33,dof21) + ke(9,4); % Col 4 - element dof21
    K(dof33,dof22) = K(dof33,dof22) + ke(9,5); % Col 5 - element dof22
    K(dof33,dof23) = K(dof33,dof23) + ke(9,6); % Col 6 - element dof23
    K(dof33,dof31) = K(dof33,dof31) + ke(9,7); % Col 7 - element dof31
    K(dof33,dof32) = K(dof33,dof32) + ke(9,8); % Col 8 - element dof32
    K(dof33,dof33) = K(dof33,dof33) + ke(9,9); % Col 9 - element dof33
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SOLVER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specification of submatrices directly
KRR = K(dofs_restrained,dofs_restrained);
KRF = K(dofs_restrained,dofs_free);
KFR = K(dofs_free,dofs_restrained);
KFF = K(dofs_free,dofs_free);
fF = F(dofs_free);
% Cannot yet form fR as the reactions are unknown

% Solution for the unknown nodal dofs
uR = zeros(length(dofs_restrained),1); % BC - zero displacement on restricted nodes
uF = KFF\(fF - KFR*uR); % 1st matrix equation
U = zeros(3*nodes,1);
U(dofs_restrained) = uR;
U(dofs_free) = uF; % full nodal dof vector

% Solution for the unknown reactions
fR = KRF*uF + KRR*uR; % 2nd matrix equation
F(dofs_restrained) = fR; % full nodal force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% POST-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating nodal coordinates after deformation
% Creating two vectors - one with real size deformed coordinates and one
% with amplified deformed coordinate (for plotting purposes)
NODES.new_coords = zeros(size(NODES.coords)); 
NODES.amp_coords = zeros(size(NODES.coords)); 
amp = 1; % amplification factor for plotting purposes only 
for I = 1:size(NODES.coords,1)
    for J = 1:size(NODES.coords,2)    
        NODES.amp_coords(I,J) = NODES.coords(I,J) + U(NODES.dofs(I,J))*amp;
        NODES.new_coords(I,J) = NODES.coords(I,J) + U(NODES.dofs(I,J));
    end
end

% Plotting
figure('units','normalized','outerposition',[0 0 1 1]); hold all; grid on; tol = 1e-3;
xmin = min(NODES.amp_coords(:,1)); xmax = max(NODES.amp_coords(:,1)); difx = xmax - xmin;
ymin = min(NODES.amp_coords(:,2)); ymax = max(NODES.amp_coords(:,2)); dify = ymax - ymin; fac = 0.25;
% axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]);
% Note that if the 'squished' structural shape bothers you, replace the
% above line with 'axis equal'
axis equal;

for EL = 1:elements
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); % identify element node numbers
    
%     Plotting original structure
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y original coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y original coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y original coordinates

    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
    
    patch([x1,x2,x3],[y1,y2,y3],[0.5 0.5 0.5]);
    
    
%     Check on changes in member lengths and plotting amplified deformed structure
    x1_amp = NODES.amp_coords(n1,1); y1_amp = NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
    x2_amp = NODES.amp_coords(n2,1); y2_amp = NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
    x3_amp = NODES.amp_coords(n3,1); y3_amp = NODES.amp_coords(n3,2); % element node 3 - x,y amplified deformed coordinates
    
    patch([x1_amp,x2_amp,x3_amp],[y1_amp,y2_amp,y3_amp],'r');
    
    
%     Plotting nodes last!
%     plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x3,y3,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');
%     plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
%     plot(x3_amp,y3_amp,'ko','Markersize',7,'MarkerFaceColor','y');
end
xlabel('x coordinate');
ylabel('y coordinate');
title('Deformed shape')
set(gca,'FontSize',30);