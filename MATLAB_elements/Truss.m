% Iñigo Basterretxea Jacob 01246662

% Script that implements the Finite Element Method for a rectangular domain
% with truss elements
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSOR MODULE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
E = 2e5; % N/mm2 - modulus of elasticity of each bar element
A = 100*pi; % mm2 - cross-sectional area of each bar element
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

elements = nx*3*ny + nx + ny;
dx = Lx/nx; % Distance between nodes in x
dy = Ly/ny; % Distance between nodes in y
Nx = nx+1; % Nodes in x direction
Ny = ny+1; % Nodes in y direction
nodes = Nx*Ny; % Total number of nodes

% Note that node numbers are identified by their row number
% Specifying nodal x-y coordinates. 1st column is x & 2nd is y
for i=1:Nx
    for j=1:Ny
        NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
    end
end

% Specifying  DOFs in pairs per node. Odd numbers are horizontal DOFs &
% even numbers are vertical DOFs
NODES.dofs = zeros(nodes,2);
for i=1:nodes  
    NODES.dofs(i,:) = [2*i-1, 2*i];
end
     
% Specifying element nodal connectivity in rows of 2 nodes per element
% (order of the elements does not matter)
ELEMENTS = zeros(elements,2);
for i = 0:ny-1
    % Horizontal member
    for k = 1:nx
        ELEMENTS((3*nx+1)*i+k,:) = [i+1, i+ny+2] + [(k-1)*(ny+1), (k-1)*(ny+1)]; 
    end
    % Vertical member
    for j = 1:nx+1
        ELEMENTS((3*nx+1)*i+nx+j,:) = [i+1, i+2] + [(j-1)*(ny+1), (j-1)*(ny+1)]; 
    end
    % High-low diagonal
    for l = 1:nx
        ELEMENTS((3*nx+1)*i+2*nx+1+l,:) = [i+2, i+ny+2] + [(l-1)*(ny+1), (l-1)*(ny+1)]; 
    end
end
% Top horizontal members
for k = 1:nx
        ELEMENTS((3*nx+1)*ny+k,:) = [ny+1, ny+ny+2] + [(k-1)*(ny+1), (k-1)*(ny+1)];
end

% Degrees of freedom & other parameters
% unknown nodal x-y dofs
dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*Ny-2,2*Ny:2*nodes]; 
% known nodal x-y dofs due to BC
dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*Ny-1];

count = 1; % Initialise q-vector counter
q = loadw.*ones(nx,1); % Vector of upper boundary UDLs
F = zeros(2*nodes,1); % Initialise nodal forces vector

K = zeros(2*nodes); % Initialise global stiffness matrix

% Constructing the global stiffness matrix and the nodal forces vector
for EL = 1:elements % loop through all elements & build stiffness matrix
    
    % identify element node numbers
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2);
    
    % identify node coordinates
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2);
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2);
    
    % identify element DOFs
    dof11 = NODES.dofs(n1,1); dof12 = NODES.dofs(n1,2);
    dof21 = NODES.dofs(n2,1); dof22 = NODES.dofs(n2,2);
    x21 = x2 - x1;
    
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE F %%%
    %%%%%%%%%%%%%%%%%%
    % Element length
    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); 
    
    % Weight contribution at all elements
    F(dof12) = F(dof12) + A*weight*Le*1e-9/3;
    F(dof22) = F(dof22) + A*weight*Le*1e-9/3;
    
    % UDL contribution at top boundary elements
    if y1 == Ly  && y2 == Ly
        F(dof12) = F(dof12) + q(count)*abs(x21)*1e-3/2;
        F(dof22) = F(dof22) + q(count)*abs(x21)*1e-3/2;
        count = count + 1;
    end
    
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE K %%%
    %%%%%%%%%%%%%%%%%%
    
    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE direction of the x axis
    c = cos(alpha); c2 = c*c; s = sin(alpha); s2 = s*s; cs = c*s; % angle parameters
    ke = E*A/Le; % element axial stiffness
    
    % Updating global stiffness matrix [K] coefficients 
    % Note that for each element you still have to 'think locally'
    % Row 1 - element dof11
    K(dof11,dof11) = K(dof11,dof11) + ke*c2; % Col 1 - element dof11
    K(dof11,dof12) = K(dof11,dof12) + ke*cs; % Col 2 - element dof12
    K(dof11,dof21) = K(dof11,dof21) - ke*c2; % Col 3 - element dof21
    K(dof11,dof22) = K(dof11,dof22) - ke*cs; % Col 4 - element dof22
    
    % Row 2 - element dof12
    K(dof12,dof11) = K(dof12,dof11) + ke*cs; % Col 1 - element dof11
    K(dof12,dof12) = K(dof12,dof12) + ke*s2; % Col 2 - element dof12
    K(dof12,dof21) = K(dof12,dof21) - ke*cs; % Col 3 - element dof21
    K(dof12,dof22) = K(dof12,dof22) - ke*s2; % Col 4 - element dof22
    
    % Row 3 - element dof21
    K(dof21,dof11) = K(dof21,dof11) - ke*c2; % Col 1 - element dof11
    K(dof21,dof12) = K(dof21,dof12) - ke*cs; % Col 2 - element dof12
    K(dof21,dof21) = K(dof21,dof21) + ke*c2; % Col 3 - element dof21
    K(dof21,dof22) = K(dof21,dof22) + ke*cs; % Col 4 - element dof22
    
    % Row 4 - element dof22
    K(dof22,dof11) = K(dof22,dof11) - ke*cs; % Col 1 - element dof11
    K(dof22,dof12) = K(dof22,dof12) - ke*s2; % Col 2 - element dof12
    K(dof22,dof21) = K(dof22,dof21) + ke*cs; % Col 3 - element dof21
    K(dof22,dof22) = K(dof22,dof22) + ke*s2; % Col 4 - element dof22
    
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
U = zeros(2*nodes,1);
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
axis equal;
% axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]);
% Note that if the 'squished' structural shape bothers you, replace the
% above line with 'axis equal'

for EL = 1:elements
    n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); % identify element node numbers
    
    % Plotting original structure
    x1 = NODES.coords(n1,1); y1 = NODES.coords(n1,2); % element node 1 - x,y original coordinates
    x2 = NODES.coords(n2,1); y2 = NODES.coords(n2,2); % element node 2 - x,y original coordinates
    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
    ke = E*A/Le; % element axial stiffness
    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
    plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 
    
    % Check on changes in member lengths and plotting amplified deformed structure
    x1_amp = NODES.amp_coords(n1,1); y1_amp = NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
    x2_amp = NODES.amp_coords(n2,1); y2_amp = NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates     
    x1_new = NODES.new_coords(n1,1); y1_new = NODES.new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
    x2_new = NODES.new_coords(n2,1); y2_new = NODES.new_coords(n2,2); % element node 2 - x,y actual deformed coordinates    
    u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % reconstruction of element global dofs
    up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup = up2 - up1; % reconstruction of element local dofs     
    % note that this now gives you access to the element forces       
    if dup < -tol % element length has decreased - member in compression
        col = 'b'; % blue colour
    elseif dup > tol % element length as increased - member in tension
        col = 'r'; % red colour
    else % no change in element length
        col = 'k'; % black colour
    end
    plot([x1_amp,x2_amp],[y1_amp,y2_amp],col,'Linewidth',3);
    
%     % Calculate the axial force in the element, and display it within a text box
%     Fax = ke*dup;
%     x_lims = get(gca,'xlim'); xmin = x_lims(1); xmax = x_lims(2); 
%     y_lims = get(gca,'ylim'); ymin = y_lims(1); ymax = y_lims(2);
%     ax_pos = get(gca,'position'); ax_xpos = ax_pos(1); ax_ypos = ax_pos(2); ax_width = ax_pos(3); ax_height = ax_pos(4);
%     txt_x_pos = (0.5*(x1_amp + x2_amp) - xmin)/(xmax - xmin) * ax_width + ax_xpos;
%     txt_y_pos = (0.5*(y1_amp + y2_amp) - ymin)/(ymax - ymin) * ax_height + ax_ypos;
%     txt_pos = [txt_x_pos txt_y_pos 0 0];
%     annotation('textbox',txt_pos,'String',num2str(Fax),'FitBoxToText','on','color',col,'FontSize',20,'BackgroundColor','w');
    
    % Plotting nodes last!
%     plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
%     plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
%     plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
%     plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
end
xlabel('x coordinate');
ylabel('y coordinate');
set(gca,'FontSize',30);

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