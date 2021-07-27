% Using Lecture_4_2Dtruss by Dr. Sadowski as a template
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimensions of domain and elements in each direction
nx = 80;
ny = 4;
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
loadw = -1e4; % kN/m2 - Imposed UDL
q = loadw.*ones(nx,1);
q(2:2:end) = -0.5;
weight = -80e3; % kN/m3 - Unit weight of steel

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
K = zeros(2*nodes); % initialising en empty 12x12 matrix; 6 nodes at 2 dofs/node
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

% figure;
% spy(K_sparse);
% title('Sparse');
% 
% figure;
% spy(K);
% title('K');

% Get nodal force vector
forces = K_sparse*U;

% % Fitness function
% f = @(q)CST_UDL_Fitness_Single(q,t,weight,Ly,coords,elem,dofs_free,forces);
% [UDL,fval,flag,output] = fmincon(f,0,[],[]);

% disp(['Imposed UDL is = ',num2str(UDL/1e3),'kN/m']);

% Fitness function
fun = @(q) CST_UDL_Fitness(q,t,weight,Ly,coords,elem,dofs_free,forces);

x0 = zeros(nx,1);

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e6,'MaxIterations',1e4);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e10,'MaxIterations',1e5);
As = [];
bs = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];

tic;
[vars,fval2,exitflag,~,~,grad,hessian] = fmincon(fun,x0,As,bs,Aeq,beq,lb,ub,nonlcon,options);
t1 = toc;

disp(['Imposed UDL is = ',num2str(vars(1)/1e3),'N/m']);
disp(['Optimisation time = ',num2str(t1)]);

% % Note that this portion is misleading in its minimalism. We exploit the
% % fact that Matlab has efficient routines for matrix operations, and we can
% % manipulate matrices just like simple algebra. Additionally, the matrices
% % are small and the analysis is linear, so the 'inversion' of KFF is very
% % fast. Do not be fooled - for large nonlinear models this step can take
% % the lion's share of modelling time. Hours, if not days. Just ask a PhD
% % student!
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% POST-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Updating nodal coordinates after deformation
% % Creating two vectors - one with real size deformed coordinates and one
% % with amplified deformed coordinate (for plotting purposes)
% NODES.new_coords = zeros(size(NODES.coords)); 
% NODES.amp_coords = zeros(size(NODES.coords)); 
% amp = 10; % amplification factor for plotting purposes only 
% for I = 1:size(NODES.coords,1)
%     for J = 1:size(NODES.coords,2)    
%         NODES.amp_coords(I,J) = NODES.coords(I,J) + U(NODES.dofs(I,J))*amp;
%         NODES.new_coords(I,J) = NODES.coords(I,J) + U(NODES.dofs(I,J));
%     end
% end
% 
% % Plotting
% figure('units','normalized','outerposition',[0 0 1 1]); hold all; grid on; tol = 1e-3;
% xmin = min(NODES.amp_coords(:,1)); xmax = max(NODES.amp_coords(:,1)); difx = xmax - xmin;
% ymin = min(NODES.amp_coords(:,2)); ymax = max(NODES.amp_coords(:,2)); dify = ymax - ymin; fac = 0.25;
% axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]);
% % Note that if the 'squished' structural shape bothers you, replace the
% % above line with 'axis equal'
% axis equal;
% 
% for EL = 1:elements
%     n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); % identify element node numbers
%     
% %     Plotting original structure
%     x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y original coordinates
%     x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y original coordinates
%     x3 = NODES.coords(n3,1); y3 = NODES.coords(n3,2); % element node 3 - x,y original coordinates
% 
%     alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
%     patch([x1,x2,x3],[y1,y2,y3],[0.5 0.5 0.5]); 
%     
% %     Check on changes in member lengths and plotting amplified deformed structure
%     x1_amp = NODES.amp_coords(n1,1); y1_amp = NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
%     x2_amp = NODES.amp_coords(n2,1); y2_amp = NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
%     x3_amp = NODES.amp_coords(n3,1); y3_amp = NODES.amp_coords(n3,2); % element node 3 - x,y amplified deformed coordinates
% 
%     patch([x1_amp,x2_amp,x3_amp],[y1_amp,y2_amp,y3_amp],'r');
%     
%     
% %     Plotting nodes last!
% %     plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
% %     plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
% %     plot(x3,y3,'ko','Markersize',7,'MarkerFaceColor','w');
% %     plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
% %     plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
% %     plot(x3_amp,y3_amp,'ko','Markersize',7,'MarkerFaceColor','y');
% end
% xlabel('x coordinate');
% ylabel('y coordinate');
% title('Deformed shape')
% set(gca,'FontSize',30);
% 
% % Printing computed dofs & reactions - this is an important part of the post-processing,
% % make sure to include something like this in every analysis that you do.
% for dof = 1:length(uF)
%     disp(['The value of dof ',num2str(dofs_free(dof)),' is ',num2str(uF(dof))]);    
% end
% disp(' ');
% for react = 1:length(fR)
%     disp(['The value of the reaction at dof ',num2str(dofs_restrained(react)),' is ',num2str(fR(react))]);
% end
% disp(' '); disp('Vertical equilibrium check:');
% disp(['Total vertical reactions = ',num2str(fR(2))]);
% disp(['Total applied vertical loads = ',num2str(sum(F(2*Ny:2*Ny:2*N)))]);
% if abs(fR(2)+sum(F(2*Ny:2*Ny:2*N))) < 1e-6; disp('Ok.'); end
% disp(' '); disp('Horizontal equilibrium check:');
% disp(['Total horizontal reactions = ',num2str(fR(1)+fR(3))]);
% disp('Total applied horizontal loads = 0');
% if abs(fR(1)+fR(3)) < 1e-6; disp('Ok.'); end