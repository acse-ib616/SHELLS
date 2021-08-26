% Iñigo Basterretxea Jacob 01246662

% This function assembles the nodal forces vector for the truss element
% and calculates the residual R = y - f(x) with respect to the input y at the FREE DOFs

% Inputs:
% - q: vector of top boundary UDLs
% - A: cross-sectional area of element
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
function [fitness] = Truss_Fitness(q,weight,A,Ly,elem,coords,dofs,dofs_free,forces)

elements = length(elem)/2;  % no. of elements

% Initialise nodal forces vector
fuerzas = zeros(length(forces),1);
% Initialse counter for vector of top boundary UDLs
count = 1;
% Constructing the nodal forces vector
for EL = 1:elements % loop through all elements & build stiffness matrix
    
    % identify element node numbers
    n1 = elem((EL-1)*2+1); n2 = elem((EL-1)*2+2);
    
    % identify node coordinates
    x1 = coords((n1-1)*2+1); y1 = coords((n1-1)*2+2);
    x2 = coords((n2-1)*2+1); y2 = coords((n2-1)*2+2);
    
    % identify element DOFs
    dof12 = dofs(n1,2);
    dof22 = dofs(n2,2);
    
    %%%%%%%%%%%%%%%%%%
    %%% ASSEMBLE F %%%
    %%%%%%%%%%%%%%%%%%
    % Element length
    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 );
    
    % Length projected in x
    x21 = x2 - x1;
    
    % Weight contribution at all elements
    fuerzas(dof12) = fuerzas(dof12) + A*weight*Le*1e-9/3;
    fuerzas(dof22) = fuerzas(dof22) + A*weight*Le*1e-9/3;
    
    % UDL contribution at top boundary elements
    if y1 == Ly  && y2 == Ly
        fuerzas(dof12) = fuerzas(dof12) + q(count)*abs(x21)*1e-3/2;
        fuerzas(dof22) = fuerzas(dof22) + q(count)*abs(x21)*1e-3/2;
        count = count + 1;
    end
    
end


% Calculate residual ignoring restrained dofs (i.e. reactions)
fitness =  forces(dofs_free) - fuerzas(dofs_free);

end