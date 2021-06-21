% Using Lecture_4_2Dtruss by Dr. Sadowski as a template
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
E = 42e9; % N/m2 - modulus of elasticity of each bar element
nu = 0.2; % Poisson coefficient
t = 0.05; % Element thickness in m

% Specifying nodal x-y coordinates
NODES.coords = [0   0 ; % node 1 - x,y coordinates
                3   2 ; % node 2 - x,y coordinates
                6   0 ; % node 3 - x,y coordinates
                6   5 ; % node 4 - x,y coordinates
                3   5 ;
                0   5]; % node 5 - x,y coordinates
NODES.dofs = [1 2 ;  % node 1 - dofs u1,v1 (1,2) - FIXED  
              3 4 ;  % node 2 - dofs u2,v2 (3,4) - FIXED
              5 6 ;  % node 3 - dofs u3,v3 (5,6) - FIXED
              7 8 ;  % node 4 - dofs u4,v4 (7,8) - FIXED
              9 10;  % node 5 - dofs u5,v5 (9,10) - FREE
              11 12]; % node 6 - dofs u6,v6 (11,12) - FIXED
% Note that node numbers are identified by their row number
     
% Specifying element nodal connectivity (order does not matter)
ELEMENTS = [1 2 6 ; % element 1 - element dofs u1,v1,u3,v3 (1,2,5,6)
            2 3 4 ; % element 2 - element dofs u2,v2,u3,v3 (3,4,5,6)
            2 4 5 ; % element 3 - element dofs u2,v2,u4,v4 (3,4,7,8)
            2 5 6 ]; % element 4 - element dofs u3,v3,u4,v4 (5,6,7,8)

% Degrees of freedom & other parameters
dofs_free = [9,10]; % unknown nodal x-y dofs (at nodes 3,4,5)
dofs_restrained = [1,2,3,4,5,6,7,8,11,12]; % known nodal x-y dofs due to BC (at nodes 1,2)
nodes = size(NODES.coords,1); % no. of nodes i.e. 6
elements = size(ELEMENTS,1); % no. of elements i.e.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the global stiffness matrix
K = zeros(2*nodes); % initialising en empty 12x12 matrix; 6 nodes at 2 dofs/node
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); % identify element node numbers
    n3 = ELEMENTS(EL,3);
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y coordinates
    
%     x21 = x2 - x1; x31 = x3 - x1; % Triangle sides
%     y21 = y2 - y1; y31 = y3 - y1;
%     Jacobian = abs(x21*y31 - x31*y21);
%     A = Jacobian/2; % m2 - cross-sectional area of triangular element
    
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2); % element node 1 - dofs
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2); % element node 2 - dofs
    dof31 = NODES.dofs(n3,1); dof32 = NODES.dofs(n3,2); % element node 3 - dofs
    
    constants = E*t/(1-nu^2); % Plane stress condition
    
    
%     y23 = y2 - y3; y12 = -y21; x32 = x3 - x2; x13 = -x31; % Triangle sides
%     
%     B = 1/(2*A)*[y23 0 y31 0 y12 0;
%                  0 x32 0 x13 0 x21;
%                  x32 y23 x13 y31 x21 y12];

    [r,wr] = lgwt(2,0.0,1.0); % Gauss quadrature
%     [r,wr] = md_gauss(2);
    s = r; ws = wr;
    
    ke = zeros(6,6);
    
    for m=1:length(r)
        for n=1:length(s)
            kers = k_CST(x1,x2,x3,y1,y2,y3,nu,r(m),s(n));   
            
            ke = ke+kers.*wr(m).*ws(n);
        end
    end
    
    
    ke = triu(ke)+triu(ke,1)'; % mirror upper matrix to the lower half

    ke = 0.5*constants*ke; % element stiffness. The 0.5 comes from integrating the area of a triangle (bh/2)
    
    % Updating global stiffness matrix [K] coefficients 
    % Note that for each element you still have to 'think locally'
    % Row 1 - element dof11
    K(dof11,dof11) = K(dof11,dof11) + ke(1,1); % Col 1 - element dof11
    K(dof11,dof12) = K(dof11,dof12) + ke(1,2); % Col 2 - element dof12
    K(dof11,dof21) = K(dof11,dof21) + ke(1,3); % Col 3 - element dof21
    K(dof11,dof22) = K(dof11,dof22) + ke(1,4); % Col 4 - element dof22
    K(dof11,dof31) = K(dof11,dof31) + ke(1,5); % Col 5 - element dof31
    K(dof11,dof32) = K(dof11,dof32) + ke(1,6); % Col 6 - element dof32
    
    % Row 2 - element dof12
    K(dof12,dof11) = K(dof12,dof11) + ke(2,1); % Col 1 - element dof11
    K(dof12,dof12) = K(dof12,dof12) + ke(2,2); % Col 2 - element dof12
    K(dof12,dof21) = K(dof12,dof21) + ke(2,3); % Col 3 - element dof21
    K(dof12,dof22) = K(dof12,dof22) + ke(2,4); % Col 4 - element dof22
    K(dof12,dof31) = K(dof12,dof31) + ke(2,5); % Col 5 - element dof31
    K(dof12,dof32) = K(dof12,dof32) + ke(2,6); % Col 6 - element dof32
    
    % Row 3 - element dof21
    K(dof21,dof11) = K(dof21,dof11) + ke(3,1); % Col 1 - element dof11
    K(dof21,dof12) = K(dof21,dof12) + ke(3,2); % Col 2 - element dof12
    K(dof21,dof21) = K(dof21,dof21) + ke(3,3); % Col 3 - element dof21
    K(dof21,dof22) = K(dof21,dof22) + ke(3,4); % Col 4 - element dof22
    K(dof21,dof31) = K(dof21,dof31) + ke(3,5); % Col 5 - element dof31
    K(dof21,dof32) = K(dof21,dof32) + ke(3,6); % Col 6 - element dof32
    
    % Row 4 - element dof22
    K(dof22,dof11) = K(dof22,dof11) + ke(4,1); % Col 1 - element dof11
    K(dof22,dof12) = K(dof22,dof12) + ke(4,2); % Col 2 - element dof12
    K(dof22,dof21) = K(dof22,dof21) + ke(4,3); % Col 3 - element dof21
    K(dof22,dof22) = K(dof22,dof22) + ke(4,4); % Col 4 - element dof22
    K(dof22,dof31) = K(dof22,dof31) + ke(4,5); % Col 5 - element dof31
    K(dof22,dof32) = K(dof22,dof32) + ke(4,6); % Col 6 - element dof32
    
    % Row 5 - element dof31
    K(dof31,dof11) = K(dof31,dof11) + ke(5,1); % Col 1 - element dof11
    K(dof31,dof12) = K(dof31,dof12) + ke(5,2); % Col 2 - element dof12
    K(dof31,dof21) = K(dof31,dof21) + ke(5,3); % Col 3 - element dof21
    K(dof31,dof22) = K(dof31,dof22) + ke(5,4); % Col 4 - element dof22
    K(dof31,dof31) = K(dof31,dof31) + ke(5,5); % Col 5 - element dof31
    K(dof31,dof32) = K(dof31,dof32) + ke(5,6); % Col 6 - element dof32
    
    % Row 6 - element dof32
    K(dof32,dof11) = K(dof32,dof11) + ke(6,1); % Col 1 - element dof11
    K(dof32,dof12) = K(dof32,dof12) + ke(6,2); % Col 2 - element dof12
    K(dof32,dof21) = K(dof32,dof21) + ke(6,3); % Col 3 - element dof21
    K(dof32,dof22) = K(dof32,dof22) + ke(6,4); % Col 4 - element dof22
    K(dof32,dof31) = K(dof32,dof31) + ke(6,5); % Col 5 - element dof31
    K(dof32,dof32) = K(dof32,dof32) + ke(6,6); % Col 6 - element dof32
    
    
end

% Constructing the global nodal force vector
F = zeros(2*nodes,1); % initialising en empty a 12x1 column vector for convenience
F(10) = -2250e3; F(8) = -750e3; % The load P acts downwards on node 5 i.e. it affects global dof 10
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
uR = zeros(length(dofs_restrained),1); % BC - zero displacement on nodes 1,2,3,4,5,6,7,8,11,12
uF = KFF\(fF - KFR*uR); % 1st matrix equation
U = zeros(2*nodes,1);
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
amp = 200; % amplification factor for plotting purposes only 
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
axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]);
% Note that if the 'squished' structural shape bothers you, replace the
% above line with 'axis equal'

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
    plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
    plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
    plot(x3,y3,'ko','Markersize',7,'MarkerFaceColor','w');
    plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
    plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
    plot(x3_amp,y3_amp,'ko','Markersize',7,'MarkerFaceColor','y');
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
disp(['Total vertical reactions = ',num2str(F(2) + F(4) + F(6) + F(8) + F(12))]);
disp(['Total applied vertical loads = ',num2str(F(10))]);
if abs(F(2) + F(4) + F(6) + F(8) + F(10) + F(12)) < 1e-6; disp('Ok.'); end
disp(' '); disp('Horizontal equilibrium check:');
disp(['Total horizontal reactions = ',num2str(F(1) + F(3) + F(5) + F(7) + F(11))]);
disp('Total applied horizontal loads = 0');
if abs(F(1) + F(3) + F(5) + F(7) + F(9) + F(11)) < 1e-6; disp('Ok.'); end