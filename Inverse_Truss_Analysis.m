% Iñigo Basterretxea Jacob
% 01246662
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
E = 2e5; % N/mm2 - modulus of elasticity of each bar element
A = 100*pi; % mm2 - cross-sectional area of each bar element
weight = -80e3; % N/m3 - Unit weight of steel

% Dimensions of domain and elements in each direction
nx = 20*4;
ny = 1*4;
n_el = nx*3*ny + nx + ny;
Lx = 2000;
Ly = 100;
dx = Lx/nx; % Distance between nodes in x
dy = Ly/ny; % Distance between nodes in y
Nx = nx+1; % Nodes in x direction
Ny = ny+1; % Nodes in y direction
N = Nx*Ny; % Total number of nodes

loadw = -1e5; % N/m2 - Imposed UDL

q_indices = randi([0 1],nx,1);
q_values = loadw.*rand(nx,1);
q = q_values.*q_indices;
% q = loadw*ones(nx,1);

% Specifying nodal x-y coordinates
for i=1:Nx
    for j=1:Ny
        NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
    end
end

NODES.dofs = zeros(N,2);
for i=1:N  
    NODES.dofs(i,:) = [2*i-1, 2*i];
end
     
% Specifying element nodal connectivity (order does not matter)
ELEMENTS = zeros(n_el,2);
for i = 0:ny-1
    for k = 1:nx
        ELEMENTS((3*nx+1)*i+k,:) = [i+1, i+ny+2] + [(k-1)*(ny+1), (k-1)*(ny+1)]; %Horizontal member
    end
    for j = 1:nx+1
        ELEMENTS((3*nx+1)*i+nx+j,:) = [i+1, i+2] + [(j-1)*(ny+1), (j-1)*(ny+1)]; %Vertical member
    end
    for l = 1:nx
        ELEMENTS((3*nx+1)*i+2*nx+1+l,:) = [i+2, i+ny+2] + [(l-1)*(ny+1), (l-1)*(ny+1)]; %High-low diagonal
    end
end
% Top horizontal members
for k = 1:nx
        ELEMENTS((3*nx+1)*ny+k,:) = [ny+1, ny+ny+2] + [(k-1)*(ny+1), (k-1)*(ny+1)];
end

% Degrees of freedom & other parameters
dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*Ny-2,2*Ny:2*N]; % unknown nodal x-y dofs
dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*Ny-1]; % known nodal x-y dofs due to BC (at nodes 1,9)
nodes = size(NODES.coords,1); % no. of nodes i.e. 5
elements = size(ELEMENTS,1); % no. of elements i.e. 6
EA = E*A; % for convenience, solely because E and A are constant in this example


F = zeros(2*nodes,1); % initialising en empty a 10x1 column vector for convenience
count = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ASSEMBLER MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the nodal equivalent loads vector
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); % identify element node numbers
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y coordinates
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2); % element node 1 - dofs
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2); % element node 2 - dofs
    x21 = x2 - x1;
    
    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
    
    F(dof12) = F(dof12) + A*weight*Le*1e-9/3; % Weight contribution
    F(dof22) = F(dof22) + A*weight*Le*1e-9/3;
    
    if y1 == Ly  && y2 == Ly % UDL contribution
        F(dof12) = F(dof12) + q(count)*abs(x21)*1e-3/2;
        F(dof22) = F(dof22) + q(count)*abs(x21)*1e-3/2;
        count = count + 1;
    end
    
end

coords = reshape(NODES.coords',1,[]);
elem = uint32(reshape(ELEMENTS',1,[]));


[k_values,k_rows,k_cols] = Truss_K(E,A,coords,elem);
k_rows = double(k_rows) + ones(size(k_rows));
k_cols = double(k_cols) + ones(size(k_cols));
K = sparse(k_rows,k_cols,k_values,2*nodes,2*nodes);

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
uR = zeros(length(dofs_restrained),1); % BC - zero displacement on nodes 1 & 2
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

forces = K*U;

x0 = zeros(length(q),1);

tic;
[vars,fval,rel_tol,it] = GN_Truss(x0,Ly,weight,A,elem,coords,NODES.dofs,dofs_free,F,1e-4,1e2,0.0);
t1 = toc;

disp(['Imposed Load is = ',num2str(vars(1)/1e3),'kN']);
disp(['Optimisation time = ',num2str(t1)]);


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
% 
% % Plotting
% figure('units','normalized','outerposition',[0 0 1 1]); hold all; grid on; tol = 1e-3;
% xmin = min(NODES.amp_coords(:,1)); xmax = max(NODES.amp_coords(:,1)); difx = xmax - xmin;
% ymin = min(NODES.amp_coords(:,2)); ymax = max(NODES.amp_coords(:,2)); dify = ymax - ymin; fac = 0.25;
% axis equal;
% % axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]);
% % Note that if the 'squished' structural shape bothers you, replace the
% % above line with 'axis equal'
% 
% for EL = 1:elements
%     n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); % identify element node numbers
%     
%     % Plotting original structure
%     x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y original coordinates
%     x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y original coordinates
%     Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
%     ke = EA/Le; % element axial stiffness
%     alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
%     plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 
%     
%     % Check on changes in member lengths and plotting amplified deformed structure
%     x1_amp = NODES.amp_coords(n1,1); y1_amp = NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
%     x2_amp = NODES.amp_coords(n2,1); y2_amp = NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates     
%     x1_new = NODES.new_coords(n1,1); y1_new = NODES.new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
%     x2_new = NODES.new_coords(n2,1); y2_new = NODES.new_coords(n2,2); % element node 2 - x,y actual deformed coordinates    
%     u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % reconstruction of element global dofs
%     up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup = up2 - up1; % reconstruction of element local dofs     
%     % note that this now gives you access to the element forces       
% %     if dup < -tol % element length has decreased - member in compression
% %         col = 'b'; % blue colour
% %     elseif dup > tol % element length as increased - member in tension
% %         col = 'r'; % red colour
% %     else % no change in element length
% %         col = 'k'; % black colour
% %     end
% %     plot([x1_amp,x2_amp],[y1_amp,y2_amp],col,'Linewidth',3);
%     
% %     % Calculate the axial force in the element, and display it within a text box
% %     Fax = ke*dup;
% %     x_lims = get(gca,'xlim'); xmin = x_lims(1); xmax = x_lims(2); 
% %     y_lims = get(gca,'ylim'); ymin = y_lims(1); ymax = y_lims(2);
% %     ax_pos = get(gca,'position'); ax_xpos = ax_pos(1); ax_ypos = ax_pos(2); ax_width = ax_pos(3); ax_height = ax_pos(4);
% %     txt_x_pos = (0.5*(x1_amp + x2_amp) - xmin)/(xmax - xmin) * ax_width + ax_xpos;
% %     txt_y_pos = (0.5*(y1_amp + y2_amp) - ymin)/(ymax - ymin) * ax_height + ax_ypos;
% %     txt_pos = [txt_x_pos txt_y_pos 0 0];
% %     annotation('textbox',txt_pos,'String',num2str(Fax),'FitBoxToText','on','color',col,'FontSize',20,'BackgroundColor','w');
%     
%     % Plotting nodes last!
% %     plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
% %     plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
% %     plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
% %     plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
% end
% xlabel('x coordinate');
% ylabel('y coordinate');
% set(gca,'FontSize',30);

% Printing computed dofs & reactions - this is an important part of the post-processing,
% make sure to include something like this in every analysis that you do.
% for dof = 1:length(uF)
%     disp(['The value of dof ',num2str(dofs_free(dof)),' is ',num2str(uF(dof))]);    
% end
% disp(' ');
% for react = 1:length(fR)
%     disp(['The value of the reaction at dof ',num2str(dofs_restrained(react)),' is ',num2str(fR(react))]);
% end
% disp(' '); disp('Vertical equilibrium check:');
% disp(['Total vertical reactions = ',num2str(fR(2) + fR(4))]);
% disp(['Total applied vertical loads = ',num2str(sum(Point_Loads))]);
% if abs(fR(2) + fR(4) + sum(Point_Loads)) < 1e-6; disp('Ok.'); end
% disp(' '); disp('Horizontal equilibrium check:');
% disp(['Total horizontal reactions = ',num2str(fR(1) + fR(3))]);
% disp('Total applied horizontal loads = 0');
% if abs(fR(1) + fR(3)) < 1e-6; disp('Ok.'); end