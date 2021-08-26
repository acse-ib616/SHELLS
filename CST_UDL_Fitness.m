% Iñigo Basterretxea Jacob 01246662

% This function assembles the nodal forces vector for the CST element
% and calculates the residual R = y - f(x) with respect to the input F at the FREE DOFs

% Inputs:
% - q: vector of top boundary UDLs
% - t: thickness of element
% - weight: unit weight of material
% - Ly: height of domain
% - coords: vector of x,y nodal coordinates
% - elem: vector of element nodal connectivity
% - dofs: vector of DOFs of system
% - dofs_free: vector of unrestricted DOFs of system
% - forces: input nodal forces vector F

% Output:
% - fitness: vector containing residual of assembled F with respect to
% input F
function [fitness] = CST_UDL_Fitness(q,t,weight,Ly,coords,elem,dofs,dofs_free,forces)

elements = length(elem)/3;  % no. of elements

% Initialise nodal forces vector
fuerzas = zeros(length(forces),1);
% Initialse counter for vector of top boundary UDLs
count = 1;
% Constructing the nodal forces vector
for EL = 1:elements
    
    % identify element node numbers
    n1 = elem((EL-1)*3+1); n2 = elem((EL-1)*3+2); n3 = elem((EL-1)*3+3);
    
    % identify node coordinates
    x1 = coords((n1-1)*2+1); y1 = coords((n1-1)*2+2); 
    x2 = coords((n2-1)*2+1); y2 = coords((n2-1)*2+2);
    x3 = coords((n3-1)*2+1); y3 = coords((n3-1)*2+2);
    
    % identify element DOFs
    dof12 = dofs(n1,2);
    dof22 = dofs(n2,2);
    dof32 = dofs(n3,2);
    
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE F %%%
    %%%%%%%%%%%%%%%%%%
    % Triangle sides
    x21 = x2 - x1; x31 = x3 - x1; 
    y21 = y2 - y1; y31 = y3 - y1;
    
    % Area of element
    A = abs(x21*y31 - x31*y21)/2;
    
    % Weight contribution at all elements
    fuerzas(dof12) = fuerzas(dof12) + A*weight*t*1e-9/3; % Weight contribution
    fuerzas(dof22) = fuerzas(dof22) + A*weight*t*1e-9/3;
    fuerzas(dof32) = fuerzas(dof32) + A*weight*t*1e-9/3;
    
    % UDL contribution at top boundary elements
    if y1 == Ly  && y3 == Ly
        fuerzas(dof12) = fuerzas(dof12) + q(count)*abs(x31)*1e-3/2;
        fuerzas(dof32) = fuerzas(dof32) + q(count)*abs(x31)*1e-3/2;
        count = count + 1;
    end
    
end

% Calculate remainder ignoring residual dofs (i.e. reactions)
fitness = forces(dofs_free) - fuerzas(dofs_free);
end