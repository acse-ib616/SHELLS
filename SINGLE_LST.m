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
t = 0.005; % Element thickness in m


% Specifying nodal x-y coordinates
NODES.coords = [0 0;
                1 0;
                0 1;
                0.5 0;
                0.5 0.5;
                0 0.5];

NODES.dofs = [1 2;
              3 4;
              5 6;
              7 8;
              9 10;
              11 12];
% Note that node numbers are identified by their row number
     
% Specifying element nodal connectivity (order does not matter)
ELEMENTS = [1 2 3 4 5 6];

% Degrees of freedom & other parameters
dofs_free = [3,5:12]; % unknown nodal x-y dofs
dofs_restrained = [1,2,4]; % known nodal x-y dofs due to BC (at nodes 1,9)
nodes = size(NODES.coords,1); % no. of nodes i.e. 52
elements = size(ELEMENTS,1); % no. of elements i.e. 16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the global stiffness matrix
K = zeros(2*nodes); % initialising en empty 66x66 matrix; 33 nodes at 2 dofs/node
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
    
    
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2);  % element node 1 - dofs
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2);  % element node 2 - dofs
    dof31 = NODES.dofs(n3,1); dof32 = NODES.dofs(n3,2);  % element node 3 - dofs
    dof41 = NODES.dofs(n4,1); dof42 = NODES.dofs(n4,2);  % element node 4 - dofs
    dof51 = NODES.dofs(n5,1); dof52 = NODES.dofs(n5,2);  % element node 5 - dofs
    dof61 = NODES.dofs(n6,1); dof62 = NODES.dofs(n6,2);  % element node 6 - dofs
    
    [r,wr] = lgwt(3,0.0,1.0); % Gauss quadrature
%     [r,wr] = md_gauss(2);
    s = r; ws = wr;
        
    constants = E*t/(1-nu.^2); % Constants over element
    
    ke = k_mat(x1,x2,x3,y1,y2,y3,nu);
%     ke = zeros(12,12);
%     
%     for m=1:length(r)
%         for n=1:length(s)
%             kers = k_mat2(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,nu,r(m),s(n));   
%             
%             ke = ke+kers.*wr(m).*ws(n);              
%         end
%     end
    
    
%     ke = triu(ke)+triu(ke,1)'; % mirror upper matrix to the lower half

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
    K(dof11,dof41) = K(dof11,dof41) + ke(1,7); % Col 7 - element dof41
    K(dof11,dof42) = K(dof11,dof42) + ke(1,8); % Col 8 - element dof42
    K(dof11,dof51) = K(dof11,dof51) + ke(1,9); % Col 9 - element dof51
    K(dof11,dof52) = K(dof11,dof52) + ke(1,10); % Col 10 - element dof52
    K(dof11,dof61) = K(dof11,dof61) + ke(1,11); % Col 11 - element dof61
    K(dof11,dof62) = K(dof11,dof62) + ke(1,12); % Col 12 - element dof62
    
    % Row 2 - element dof12
    K(dof12,dof11) = K(dof12,dof11) + ke(2,1); % Col 1 - element dof11
    K(dof12,dof12) = K(dof12,dof12) + ke(2,2); % Col 2 - element dof12
    K(dof12,dof21) = K(dof12,dof21) + ke(2,3); % Col 3 - element dof21
    K(dof12,dof22) = K(dof12,dof22) + ke(2,4); % Col 4 - element dof22
    K(dof12,dof31) = K(dof12,dof31) + ke(2,5); % Col 5 - element dof31
    K(dof12,dof32) = K(dof12,dof32) + ke(2,6); % Col 6 - element dof32
    K(dof12,dof41) = K(dof12,dof41) + ke(2,7); % Col 7 - element dof41
    K(dof12,dof42) = K(dof12,dof42) + ke(2,8); % Col 8 - element dof42
    K(dof12,dof51) = K(dof12,dof51) + ke(2,9); % Col 9 - element dof51
    K(dof12,dof52) = K(dof12,dof52) + ke(2,10); % Col 10 - element dof52
    K(dof12,dof61) = K(dof12,dof61) + ke(2,11); % Col 11 - element dof61
    K(dof12,dof62) = K(dof12,dof62) + ke(2,12); % Col 12 - element dof62
    
    % Row 3 - element dof21
    K(dof21,dof11) = K(dof21,dof11) + ke(3,1); % Col 1 - element dof11
    K(dof21,dof12) = K(dof21,dof12) + ke(3,2); % Col 2 - element dof12
    K(dof21,dof21) = K(dof21,dof21) + ke(3,3); % Col 3 - element dof21
    K(dof21,dof22) = K(dof21,dof22) + ke(3,4); % Col 4 - element dof22
    K(dof21,dof31) = K(dof21,dof31) + ke(3,5); % Col 5 - element dof31
    K(dof21,dof32) = K(dof21,dof32) + ke(3,6); % Col 6 - element dof32
    K(dof21,dof41) = K(dof21,dof41) + ke(3,7); % Col 7 - element dof41
    K(dof21,dof42) = K(dof21,dof42) + ke(3,8); % Col 8 - element dof42
    K(dof21,dof51) = K(dof21,dof51) + ke(3,9); % Col 9 - element dof51
    K(dof21,dof52) = K(dof21,dof52) + ke(3,10); % Col 10 - element dof52
    K(dof21,dof61) = K(dof21,dof61) + ke(3,11); % Col 11 - element dof61
    K(dof21,dof62) = K(dof21,dof62) + ke(3,12); % Col 12 - element dof62
    
    % Row 4 - element dof22
    K(dof22,dof11) = K(dof22,dof11) + ke(4,1); % Col 1 - element dof11
    K(dof22,dof12) = K(dof22,dof12) + ke(4,2); % Col 2 - element dof12
    K(dof22,dof21) = K(dof22,dof21) + ke(4,3); % Col 3 - element dof21
    K(dof22,dof22) = K(dof22,dof22) + ke(4,4); % Col 4 - element dof22
    K(dof22,dof31) = K(dof22,dof31) + ke(4,5); % Col 5 - element dof31
    K(dof22,dof32) = K(dof22,dof32) + ke(4,6); % Col 6 - element dof32
    K(dof22,dof41) = K(dof22,dof41) + ke(4,7); % Col 7 - element dof41
    K(dof22,dof42) = K(dof22,dof42) + ke(4,8); % Col 8 - element dof42
    K(dof22,dof51) = K(dof22,dof51) + ke(4,9); % Col 9 - element dof51
    K(dof22,dof52) = K(dof22,dof52) + ke(4,10); % Col 10 - element dof52
    K(dof22,dof61) = K(dof22,dof61) + ke(4,11); % Col 11 - element dof61
    K(dof22,dof62) = K(dof22,dof62) + ke(4,12); % Col 12 - element dof62
    
    % Row 5 - element dof31
    K(dof31,dof11) = K(dof31,dof11) + ke(5,1); % Col 1 - element dof11
    K(dof31,dof12) = K(dof31,dof12) + ke(5,2); % Col 2 - element dof12
    K(dof31,dof21) = K(dof31,dof21) + ke(5,3); % Col 3 - element dof21
    K(dof31,dof22) = K(dof31,dof22) + ke(5,4); % Col 4 - element dof22
    K(dof31,dof31) = K(dof31,dof31) + ke(5,5); % Col 5 - element dof31
    K(dof31,dof32) = K(dof31,dof32) + ke(5,6); % Col 6 - element dof32
    K(dof31,dof41) = K(dof31,dof41) + ke(5,7); % Col 7 - element dof41
    K(dof31,dof42) = K(dof31,dof42) + ke(5,8); % Col 8 - element dof42
    K(dof31,dof51) = K(dof31,dof51) + ke(5,9); % Col 9 - element dof51
    K(dof31,dof52) = K(dof31,dof52) + ke(5,10); % Col 10 - element dof52
    K(dof31,dof61) = K(dof31,dof61) + ke(5,11); % Col 11 - element dof61
    K(dof31,dof62) = K(dof31,dof62) + ke(5,12); % Col 12 - element dof62
    
    % Row 6 - element dof32
    K(dof32,dof11) = K(dof32,dof11) + ke(6,1); % Col 1 - element dof11
    K(dof32,dof12) = K(dof32,dof12) + ke(6,2); % Col 2 - element dof12
    K(dof32,dof21) = K(dof32,dof21) + ke(6,3); % Col 3 - element dof21
    K(dof32,dof22) = K(dof32,dof22) + ke(6,4); % Col 4 - element dof22
    K(dof32,dof31) = K(dof32,dof31) + ke(6,5); % Col 5 - element dof31
    K(dof32,dof32) = K(dof32,dof32) + ke(6,6); % Col 6 - element dof32
    K(dof32,dof41) = K(dof32,dof41) + ke(6,7); % Col 7 - element dof41
    K(dof32,dof42) = K(dof32,dof42) + ke(6,8); % Col 8 - element dof42
    K(dof32,dof51) = K(dof32,dof51) + ke(6,9); % Col 9 - element dof51
    K(dof32,dof52) = K(dof32,dof52) + ke(6,10); % Col 10 - element dof52
    K(dof32,dof61) = K(dof32,dof61) + ke(6,11); % Col 11 - element dof61
    K(dof32,dof62) = K(dof32,dof62) + ke(6,12); % Col 12 - element dof62
    
    % Row 7 - element dof41
    K(dof41,dof11) = K(dof41,dof11) + ke(7,1); % Col 1 - element dof11
    K(dof41,dof12) = K(dof41,dof12) + ke(7,2); % Col 2 - element dof12
    K(dof41,dof21) = K(dof41,dof21) + ke(7,3); % Col 3 - element dof21
    K(dof41,dof22) = K(dof41,dof22) + ke(7,4); % Col 4 - element dof22
    K(dof41,dof31) = K(dof41,dof31) + ke(7,5); % Col 5 - element dof31
    K(dof41,dof32) = K(dof41,dof32) + ke(7,6); % Col 6 - element dof32
    K(dof41,dof41) = K(dof41,dof41) + ke(7,7); % Col 7 - element dof41
    K(dof41,dof42) = K(dof41,dof42) + ke(7,8); % Col 8 - element dof42
    K(dof41,dof51) = K(dof41,dof51) + ke(7,9); % Col 9 - element dof51
    K(dof41,dof52) = K(dof41,dof52) + ke(7,10); % Col 10 - element dof52
    K(dof41,dof61) = K(dof41,dof61) + ke(7,11); % Col 11 - element dof61
    K(dof41,dof62) = K(dof41,dof62) + ke(7,12); % Col 12 - element dof62
    
    % Row 8 - element dof42
    K(dof42,dof11) = K(dof42,dof11) + ke(8,1); % Col 1 - element dof11
    K(dof42,dof12) = K(dof42,dof12) + ke(8,2); % Col 2 - element dof12
    K(dof42,dof21) = K(dof42,dof21) + ke(8,3); % Col 3 - element dof21
    K(dof42,dof22) = K(dof42,dof22) + ke(8,4); % Col 4 - element dof22
    K(dof42,dof31) = K(dof42,dof31) + ke(8,5); % Col 5 - element dof31
    K(dof42,dof32) = K(dof42,dof32) + ke(8,6); % Col 6 - element dof32
    K(dof42,dof41) = K(dof42,dof41) + ke(8,7); % Col 7 - element dof41
    K(dof42,dof42) = K(dof42,dof42) + ke(8,8); % Col 8 - element dof42
    K(dof42,dof51) = K(dof42,dof51) + ke(8,9); % Col 9 - element dof51
    K(dof42,dof52) = K(dof42,dof52) + ke(8,10); % Col 10 - element dof52
    K(dof42,dof61) = K(dof42,dof61) + ke(8,11); % Col 11 - element dof61
    K(dof42,dof62) = K(dof42,dof62) + ke(8,12); % Col 12 - element dof62
    
    % Row 9 - element dof51
    K(dof51,dof11) = K(dof51,dof11) + ke(9,1); % Col 1 - element dof11
    K(dof51,dof12) = K(dof51,dof12) + ke(9,2); % Col 2 - element dof12
    K(dof51,dof21) = K(dof51,dof21) + ke(9,3); % Col 3 - element dof21
    K(dof51,dof22) = K(dof51,dof22) + ke(9,4); % Col 4 - element dof22
    K(dof51,dof31) = K(dof51,dof31) + ke(9,5); % Col 5 - element dof31
    K(dof51,dof32) = K(dof51,dof32) + ke(9,6); % Col 6 - element dof32
    K(dof51,dof41) = K(dof51,dof41) + ke(9,7); % Col 7 - element dof41
    K(dof51,dof42) = K(dof51,dof42) + ke(9,8); % Col 8 - element dof42
    K(dof51,dof51) = K(dof51,dof51) + ke(9,9); % Col 9 - element dof51
    K(dof51,dof52) = K(dof51,dof52) + ke(9,10); % Col 10 - element dof52
    K(dof51,dof61) = K(dof51,dof61) + ke(9,11); % Col 11 - element dof61
    K(dof51,dof62) = K(dof51,dof62) + ke(9,12); % Col 12 - element dof62
    
    % Row 10 - element dof52
    K(dof52,dof11) = K(dof52,dof11) + ke(10,1); % Col 1 - element dof11
    K(dof52,dof12) = K(dof52,dof12) + ke(10,2); % Col 2 - element dof12
    K(dof52,dof21) = K(dof52,dof21) + ke(10,3); % Col 3 - element dof21
    K(dof52,dof22) = K(dof52,dof22) + ke(10,4); % Col 4 - element dof22
    K(dof52,dof31) = K(dof52,dof31) + ke(10,5); % Col 5 - element dof31
    K(dof52,dof32) = K(dof52,dof32) + ke(10,6); % Col 6 - element dof32
    K(dof52,dof41) = K(dof52,dof41) + ke(10,7); % Col 7 - element dof41
    K(dof52,dof42) = K(dof52,dof42) + ke(10,8); % Col 8 - element dof42
    K(dof52,dof51) = K(dof52,dof51) + ke(10,9); % Col 9 - element dof51
    K(dof52,dof52) = K(dof52,dof52) + ke(10,10); % Col 10 - element dof52
    K(dof52,dof61) = K(dof52,dof61) + ke(10,11); % Col 11 - element dof61
    K(dof52,dof62) = K(dof52,dof62) + ke(10,12); % Col 12 - element dof62
    
    % Row 11 - element dof61
    K(dof61,dof11) = K(dof61,dof11) + ke(11,1); % Col 1 - element dof11
    K(dof61,dof12) = K(dof61,dof12) + ke(11,2); % Col 2 - element dof12
    K(dof61,dof21) = K(dof61,dof21) + ke(11,3); % Col 3 - element dof21
    K(dof61,dof22) = K(dof61,dof22) + ke(11,4); % Col 4 - element dof22
    K(dof61,dof31) = K(dof61,dof31) + ke(11,5); % Col 5 - element dof31
    K(dof61,dof32) = K(dof61,dof32) + ke(11,6); % Col 6 - element dof32
    K(dof61,dof41) = K(dof61,dof41) + ke(11,7); % Col 7 - element dof41
    K(dof61,dof42) = K(dof61,dof42) + ke(11,8); % Col 8 - element dof42
    K(dof61,dof51) = K(dof61,dof51) + ke(11,9); % Col 9 - element dof51
    K(dof61,dof52) = K(dof61,dof52) + ke(11,10); % Col 10 - element dof52
    K(dof61,dof61) = K(dof61,dof61) + ke(11,11); % Col 11 - element dof61
    K(dof61,dof62) = K(dof61,dof62) + ke(11,12); % Col 12 - element dof62
    
    % Row 12 - element dof62
    K(dof62,dof11) = K(dof62,dof11) + ke(12,1); % Col 1 - element dof11
    K(dof62,dof12) = K(dof62,dof12) + ke(12,2); % Col 2 - element dof12
    K(dof62,dof21) = K(dof62,dof21) + ke(12,3); % Col 3 - element dof21
    K(dof62,dof22) = K(dof62,dof22) + ke(12,4); % Col 4 - element dof22
    K(dof62,dof31) = K(dof62,dof31) + ke(12,5); % Col 5 - element dof31
    K(dof62,dof32) = K(dof62,dof32) + ke(12,6); % Col 6 - element dof32
    K(dof62,dof41) = K(dof62,dof41) + ke(12,7); % Col 7 - element dof41
    K(dof62,dof42) = K(dof62,dof42) + ke(12,8); % Col 8 - element dof42
    K(dof62,dof51) = K(dof62,dof51) + ke(12,9); % Col 9 - element dof51
    K(dof62,dof52) = K(dof62,dof52) + ke(12,10); % Col 10 - element dof52
    K(dof62,dof61) = K(dof62,dof61) + ke(12,11); % Col 11 - element dof61
    K(dof62,dof62) = K(dof62,dof62) + ke(12,12); % Col 12 - element dof62
    
    
end

% Constructing the global nodal force vector
F = zeros(2*nodes,1); % initialising en empty a 2Nx1 column vector for convenience
F(6) = -1e4; % The load P acts downwards on node 166 i.e. it affects global dof 166*2
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
amp = 100; % amplification factor for plotting purposes only 
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
    n4 = ELEMENTS(EL,4); n5 = ELEMENTS(EL,5); n6 = ELEMENTS(EL,6);
%     Plotting original structure
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y original coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y original coordinates
    x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y original coordinates
    x4 = NODES.coords(n4,1); y4 = NODES.coords(n4,2); % element node 4 - x,y original coordinates
    x5 = NODES.coords(n5,1); y5 = NODES.coords(n5,2); % element node 5 - x,y original coordinates
    x6 = NODES.coords(n6,1); y6 = NODES.coords(n6,2); % element node 6 - x,y original coordinates

    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
    
    patch([x1,x4,x2,x5,x3,x6],[y1,y4,y2,y5,y3,y6],[0.5 0.5 0.5]);
    
    
%     Check on changes in member lengths and plotting amplified deformed structure
    x1_amp = NODES.amp_coords(n1,1); y1_amp = NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
    x2_amp = NODES.amp_coords(n2,1); y2_amp = NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
    x3_amp = NODES.amp_coords(n3,1); y3_amp = NODES.amp_coords(n3,2); % element node 3 - x,y amplified deformed coordinates
    x4_amp = NODES.amp_coords(n4,1); y4_amp = NODES.amp_coords(n4,2); % element node 1 - x,y amplified deformed coordinates
    x5_amp = NODES.amp_coords(n5,1); y5_amp = NODES.amp_coords(n5,2); % element node 2 - x,y amplified deformed coordinates
    x6_amp = NODES.amp_coords(n6,1); y6_amp = NODES.amp_coords(n6,2); % element node 3 - x,y amplified deformed coordinates

    patch([x1_amp,x4_amp,x2_amp,x5_amp,x3_amp,x6_amp],[y1_amp,y4_amp,y2_amp,y5_amp,y3_amp,y6_amp],'r');
    
    
%     Plotting nodes last!
%     plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x4,y4,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x5,y5,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x3,y3,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x6,y6,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
%     plot(x4_amp,y4_amp,'ko','Markersize',7,'MarkerFaceColor','y');
%     plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
%     plot(x5_amp,y5_amp,'ko','Markersize',7,'MarkerFaceColor','y');
%     plot(x3_amp,y3_amp,'ko','Markersize',7,'MarkerFaceColor','y');
%     plot(x6_amp,y6_amp,'ko','Markersize',7,'MarkerFaceColor','y');
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
disp(['Total applied vertical loads = ',num2str(F(6))]);
if abs(fR(2) + fR(3) + F(6)) < 1e-6; disp('Ok.'); end
disp(' '); disp('Horizontal equilibrium check:');
disp(['Total horizontal reactions = ',num2str(fR(1))]);
disp('Total applied horizontal loads = 0');
if abs(fR(1)) < 1e-6; disp('Ok.'); end


K_info =   [1,1, 1,1,  4.375000000000001e+07;
            1,2, 1,1,  2.388135985261405e-08;
            2,1, 1,1,  1.458333333333334e+07;
            2,2, 1,1,  1.458333333333334e+07;
            3,1, 1,1, -8.769990906327194e-09;
            3,2, 1,1, -1.458333333333333e+07;
            4,1, 1,1, -5.833333333333340e+07;
            4,2, 1,1, -5.833333333333338e+07;
            5,1, 1,1, -1.303851604461670e-08;
            5,2, 1,1, -1.490116119384766e-08;
            6,1, 1,1,  8.381903171539307e-08;
            6,2, 1,1,  5.833333333333335e+07;
            1,2, 1,2,  1.093750000000000e+08;
            2,1, 1,2,  7.291666666666669e+06;
            2,2, 1,2,  3.645833333333333e+07;
            3,1, 1,2, -7.291666666666664e+06;
            3,2, 1,2, -3.507996362530876e-09;
            4,1, 1,2, -2.916666666666671e+07;
            4,2, 1,2, -1.458333333333334e+08;
            5,1, 1,2, -1.303851604461670e-08;
            5,2, 1,2, -4.097819328308105e-08;
            6,1, 1,2,  2.916666666666668e+07;
            6,2, 1,2,  6.332993507385254e-08;
            2,1, 2,1,  1.531249999999999e+08;
            2,2, 2,1,  6.562499999999999e+07;
            3,1, 2,1,  3.645833333333333e+07;
            3,2, 2,1,  1.458333333333333e+07;
            4,1, 2,1, -5.833333333333333e+07;
            4,2, 2,1, -2.916666666666665e+07;
            5,1, 2,1, -1.458333333333333e+08;
            5,2, 2,1, -5.833333333333333e+07;
            6,1, 2,1,  3.166496753692627e-08;
            6,2, 2,1,  9.778887033462524e-09;
            2,2, 2,2,  1.531249999999999e+08;
            3,1, 2,2,  7.291666666666666e+06;
            3,2, 2,2,  1.458333333333333e+07;
            4,1, 2,2, -5.833333333333331e+07;
            4,2, 2,2, -1.458333333333333e+08;
            5,1, 2,2, -2.916666666666667e+07;
            5,2, 2,2, -5.833333333333333e+07;
            6,1, 2,2,  9.313225746154785e-09;
            6,2, 2,2,  1.443549990653992e-08;
            3,1, 3,1,  1.093750000000000e+08;
            4,1, 3,1,  7.450580596923828e-09;
            4,2, 3,1, -2.602949756970297e-09;
            5,1, 3,1, -1.458333333333333e+08;
            5,2, 3,1, -2.916666666666666e+07;
            6,2, 3,1,  2.916666666666666e+07;
            3,2, 3,2,  4.374999999999999e+07;
            4,1, 3,2, -5.205899513940594e-09;
            4,2, 3,2,  1.862645149230957e-09;
            5,1, 3,2, -5.833333333333332e+07;
            5,2, 3,2, -5.833333333333333e+07;
            6,1, 3,2,  5.833333333333332e+07;
            6,2, 3,2,  3.725290298461914e-09;
            4,1, 4,1,  4.083333333333336e+08;
            4,2, 4,1,  8.750000000000006e+07;
            5,1, 4,1,  4.707245832630886e-08;
            5,2, 4,1,  8.750000000000001e+07;
            6,1, 4,1, -2.916666666666669e+08;
            6,2, 4,1, -8.750000000000004e+07;
            4,2, 4,2,  4.083333333333334e+08;
            5,1, 4,2,  8.750000000000003e+07;
            5,2, 4,2,  7.078051567077637e-08;
            6,1, 4,2, -8.750000000000003e+07;
            6,2, 4,2, -1.166666666666668e+08;
            5,1, 5,1,  4.083333333333333e+08;
            5,2, 5,1,  8.749999999999999e+07;
            6,1, 5,1, -1.166666666666667e+08;
            6,2, 5,1, -8.750000000000000e+07;
            5,2, 5,2,  4.083333333333333e+08;
            6,1, 5,2, -8.750000000000000e+07;
            6,2, 5,2, -2.916666666666667e+08;
            6,1, 6,1,  4.083333333333334e+08;
            6,2, 6,1,  8.750000000000000e+07;
            6,2, 6,2,  4.083333333333333e+08];
K2 = zeros(2*nodes);
for i=1:length(K_info)
    K2(NODES.dofs(K_info(i,1),K_info(i,2)),NODES.dofs(K_info(i,3),K_info(i,4))) = K_info(i,5);
end

U2 = [-2.85714e-4 -1e-32 0 -0.000142857 -6.42787e-20 -6.52215e-20];
U1 = [-2.85714e-4 0 5.54266e-20 -7.14286e-20 1.42857e-5 -5.71429e-5];
u = [U1(1) U2(1) U1(2) U2(2) U1(3) U2(3) U1(4) U2(4) U1(5) U2(5) U1(6) U2(6)]';