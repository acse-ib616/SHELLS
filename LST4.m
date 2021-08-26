% Iñigo Basterretxea Jacob 01246662

% Script that implements the Finite Element Method for a rectangular domain
% with LST elements
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
ny = 2; % no. of divisions along y-direction
Lx = 2000; % Height of domain in mm
Ly = 100; % Width of domain in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elements = nx*2*ny; % no. of elements
dx = Lx/(2*nx); % Distance between nodes in x
dy = Ly/(2*ny); % Distance between nodes in y
Nx = 2*nx+1; % no. of nodes in x direction
Ny = 2*ny+1; % no. of nodes in y direction
nodes = Nx*Ny; % no. of nodes

% Note that node numbers are identified by their row number
% Specifying nodal x-y coordinates. 1st column is x & 2nd is y
NODES.coords = zeros(nodes,2);
for i=1:Ny
    for j=1:Nx
        NODES.coords((i-1)*Nx+j,:) = [(j-1)*dx,dy*(i-1)];
    end
end

% Specifying  DOFs in pairs per node. Odd numbers are horizontal DOFs &
% even numbers are vertical DOFs
NODES.dofs = zeros(nodes,2);
for i=1:nodes
    
    NODES.dofs(i,:) = [2*i-1, 2*i];
end
% Note that node numbers are identified by their row number
     
% Specifying element nodal connectivity in rows of 3 nodes per element
% (order of the elements does not matter)
ELEMENTS = zeros(elements,6);
for i=1:elements
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
% unknown nodal x-y dofs
dofs_free = [1:nodes-Nx,nodes-Nx+3:2*nodes-2*Nx,2*nodes-2*Nx+2:2*nodes];
% known nodal x-y dofs due to BC
dofs_restrained = [nodes-Nx+1,nodes-Nx+2,2*nodes-2*Nx+1]; 
% dofs_free = [1:nodes-Nx,nodes-Nx+3:nodes+Nx-1,nodes+Nx+1:2*nodes]; % unknown nodal x-y dofs
% dofs_restrained = [nodes-Nx+1,nodes-Nx+2,nodes+Nx]; % known nodal x-y dofs due to BC

count = 1; % Initialise q-vector counter
q = loadw.*ones(nx,1); % Vector of upper boundary UDLs
F = zeros(2*nodes,1); % Initialise nodal forces vector

% Constructing the nodal forces vector
for EL = 1:elements % loop through all elements & build stiffness matrix
    
    % identify element node numbers
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3);
    n4 = ELEMENTS(EL,4); n5 = ELEMENTS(EL,5); n6 = ELEMENTS(EL,6);
    
    % identify node coordinates
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2);
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); 
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2);  
    x4 = NODES.coords(n4,1); y4 = NODES.coords(n4,2);  
    x5 = NODES.coords(n5,1); y5 = NODES.coords(n5,2); 
    x6 = NODES.coords(n6,1); y6 = NODES.coords(n6,2); 
    
    % identify element DOFs
    dof12 = NODES.dofs(n1,2);
    dof22 = NODES.dofs(n2,2);
    dof42 = NODES.dofs(n4,2);
    dof52 = NODES.dofs(n5,2);
    dof62 = NODES.dofs(n6,2);
    
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE F %%%
    %%%%%%%%%%%%%%%%%%
    % Triangle sides
    x21 = x2 - x1; x31 = x3 - x1; 
    y21 = y2 - y1; y31 = y3 - y1;
    
    % Half of a side
    x41 = x4 - x1;
    
    % Area of element
    A = abs(x21*y31 - x31*y21)/2; 
    
    % Weight contribution at all elements. No vertex loads!
    F(dof42) = F(dof42) + A*weight*t*1e-9/3; 
    F(dof52) = F(dof52) + A*weight*t*1e-9/3;
    F(dof62) = F(dof62) + A*weight*t*1e-9/3;
    
    % UDL contribution at top boundary elements
    if y1 == Ly  && y2 == Ly && y4 == Ly
        F(dof12) = F(dof12) + q(count)*abs(x41)*1e-3/2;
        F(dof22) = F(dof22) + q(count)*abs(x41)*1e-3/2;
        
        % Double contribution for midpoint along edge
        F(dof42) = F(dof42) + q(count)*abs(x41)*1e-3; 
        count = count + 1;
        
    end
    
end

% Make coords and elements matrices 1-dimensional. Numbers are double by
% default in MATLAB, so we must convert to unsigned integer
coords = reshape(NODES.coords',1,[]);
elem = uint32(reshape(ELEMENTS',1,[]));

% Call C++ LST stiffness matrix assembler
[k_values,k_rows,k_cols]  = LST_K(E,nu,t,coords,elem);
% Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
% position and column index vectors
k_rows = double(k_rows) + ones(size(k_rows));
k_cols = double(k_cols) + ones(size(k_cols));
K = sparse(k_rows,k_cols,k_values,2*nodes,2*nodes);


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
U = zeros(2*nodes,1);
U(dofs_restrained) = uR;
U(dofs_free) = uF; % full nodal dof vector

% Solution for the unknown reactions
fR = KRF*uF + KRR*uR; % 2nd matrix equation
F(dofs_restrained) = fR; % full nodal force vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% OPTIMISATION MODULE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess
x0 = zeros(nx,1);

% Execute Gauss-Newton algorithm
[vars,fval,rel_tol,it] = GN_LST(x0,t,weight,Ly,coords,elem,NODES.dofs,dofs_free,F,1e-4,1e6,0);

% Calculate maximum relative error
error = abs(vars - q)./abs(q).*100;
error(isinf(error)) = 0;
max_error = max(error);
if max_error < 1e-6
    disp(['CORRECT: maximum relative error is = ',num2str(max_error),'%']);
else
    disp(['ERROR: maximum relative error is larger than = ',num2str(1e-6),'%']);
end

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
axis equal;

for EL = 1:elements
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); % identify element node numbers
    
    % Plotting original structure
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y original coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y original coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y original coordinates

    patch([x1,x2,x3],[y1,y2,y3],[0.5 0.5 0.5]);
    
    
    % Check on changes in element coords and plotting amplified deformed structure
    x1_amp = NODES.amp_coords(n1,1); y1_amp = NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
    x2_amp = NODES.amp_coords(n2,1); y2_amp = NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
    x3_amp = NODES.amp_coords(n3,1); y3_amp = NODES.amp_coords(n3,2); % element node 3 - x,y amplified deformed coordinates
    
    patch([x1_amp,x2_amp,x3_amp],[y1_amp,y2_amp,y3_amp],'r');
    
end
xlabel('x coordinate');
ylabel('y coordinate');
title('Deformed shape');
set(gca,'FontSize',30);