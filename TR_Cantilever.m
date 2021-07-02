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
t = 10; % Element thickness in m

% Dimensions of domain and elements in each direction
nx = 40;
ny = 2;
n_el = nx*2*ny;
Lx = 2000;
Ly = 100;
Nx = nx+1; % Nodes in x direction (TR)
Ny = ny+1; % Nodes in y direction (TR)
N = Nx*Ny; % Total number of nodes (TR)
dx = Lx/nx; % Distance between nodes in x
dy = Ly/ny; % Distance between nodes in y

% Specifying nodal x-y coordinates
NODES.coords = zeros(N,2);
for i=1:Nx
    for j=1:Ny
        NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
    end
end

NODES.dofs = zeros(N,3);
for i=1:N
    NODES.dofs(i,:) = [3*i-2, 3*i-1, 3*i];
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


% % Degrees of freedom & other parameters
dofs_free = [1:(floor(Ny/2)+1)*3-3,(floor(Ny/2)+1)*3:N*3-(floor(Ny/2)+1)*3+1,N*3-(floor(Ny/2)+1)*3+3:3*N]; % unknown nodal x-y dofs
dofs_restrained = [(floor(Ny/2)+1)*3-2,(floor(Ny/2)+1)*3-1,N*3-(floor(Ny/2)+1)*3+2]; % known nodal x-y dofs due to BC (at nodes 1,9)
nodes = size(NODES.coords,1); % no. of nodes i.e. 52
elements = size(ELEMENTS,1); % no. of elements i.e. 16
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the global stiffness matrix
K = zeros(3*N); % initialising en empty 66x66 matrix; 33 nodes at 2 dofs/node
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); % identify element node numbers
    
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y coordinates
    
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2); dof13 = NODES.dofs(n1,3);  % element node 1 - dofs
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2); dof23 = NODES.dofs(n2,3); % element node 2 - dofs
    dof31 = NODES.dofs(n3,1); dof32 = NODES.dofs(n3,2); dof33 = NODES.dofs(n3,3); % element node 3 - dofs
    
    [r,wr] = lgwt(3,0.0,1.0); % Gauss quadrature
%     [r,wr] = md_gauss(2);
    s = r; ws = wr;
        
    constants = E*t/(1-nu.^2); % Constants over element
    
    ke_LST = zeros(12,12);
    x4 = (x1+x2)/2; y4 = (y1+y2)/2;
    x5 = (x3+x2)/2; y5 = (y3+y2)/2;
    x6 = (x3+x1)/2; y6 = (y3+y1)/2;
    
    for m=1:length(r)
        for n=1:length(s)
            kers = k_mat2(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,nu,r(m),s(n));   
            
            ke_LST = ke_LST+kers.*wr(m).*ws(n);              
        end
    end
    

    ke_LST = 0.5*constants*ke_LST; % element stiffness. The 0.5 comes from integrating the area of a triangle (bh/2)
    
    % Transformation to TR
    
    y21 = y2-y1; x21 = x2-x1; y32 = y3-y2; x32 = x3-x2;
    y13 = y1-y3; x13 = x1-x3;
    
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
    % Note that for each element you still have to 'think locally'
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

% Constructing the global nodal force vector
F = zeros(3*N,1); % initialising en empty a 2ny1 column vector for convenience
F(3*N-2) = 1e5; % The load P acts downwards on node 166 i.e. it affects global dof 166*2
F(3*N-Ny*3+1) = -1e5;
% There are no other nodal loads to apply

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
U = zeros(3*N,1);
U(dofs_restrained) = uR;
U(dofs_free) = uF; % full nodal dof vector

% Solution for the unknown reactions
fR = KRF*uF + KRR*uR; % 2nd matrix equation
F(dofs_restrained) = fR; % full nodal force vector

% Note that this portion is misleading in its minimalism. We exploit the
% fact that Matlab has efficient routines for matrix operations, and we can
% manipulate matrices just like simple algebra. Additionally, the matrices
% are small and the analysis is linear, so the 'inversion' of KFF is very
% fast. Do not be fooled - for large nonlinear models this step can take
% the lion's share of modelling time. Hours, if not days. Just ask a PhD
% student!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% POST-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating nodal coordinates after deformation
% Creating two vectors - one with real size deformed coordinates and one
% with amplified deformed coordinate (for plotting purposes)
NODES.new_coords = zeros(size(NODES.coords)); 
NODES.amp_coords = zeros(size(NODES.coords)); 
amp = 10; % amplification factor for plotting purposes only 
for I = 1:size(NODES.coords,1)
    for J = 1:size(NODES.coords,2)    
        NODES.amp_coords(I,J) = NODES.coords(I,J) + U(NODES.dofs(I,J))*amp;
        NODES.new_coords(I,J) = NODES.coords(I,J) + U(NODES.dofs(I,J));
    end
end

% Plotting
figure('units','normalized','outerposition',[0 0 1 1]); hold all; grid on; tol = 1e-3;
% xmin = min(NODES.amp_coords(:,1)); xmax = max(NODES.amp_coords(:,1)); difx = xmax - xmin;
% ymin = min(NODES.amp_coords(:,2)); ymax = max(NODES.amp_coords(:,2)); dify = ymax - ymin; fac = 0.25;
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

% Printing computed dofs & reactions - this is an important part of the post-processing,
% make sure to include something like this in every analysis that you do.
for dof = 1:length(uF)
    disp(['The value of dof ',num2str(dofs_free(dof)),' is ',num2str(uF(dof))]);    
end
disp(' ');
for react = 1:length(fR)
    disp(['The value of the reaction at dof ',num2str(dofs_restrained(react)),' is ',num2str(fR(react))]);
end
disp(' '); disp('Vertical equilibrium check:');
disp(['Total vertical reactions = ',num2str(fR(2) + fR(3))]);
disp(['Total applied vertical loads = ',num2str(0)]);
if abs(fR(2) + fR(3)) < 1e-6; disp('Ok.'); end
disp(' '); disp('Horizontal equilibrium check:');
disp(['Total horizontal reactions = ',num2str(fR(1))]);
disp(['Total applied horizontal loads = ',num2str(F(3*N-2)+F(3*N-Ny*3+1))]);
if abs(fR(1)+F(3*N-2)+F(3*N-Ny*3+1)) < 1e-6; disp('Ok.'); end