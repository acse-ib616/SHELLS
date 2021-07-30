% This function calls the MEX C++ function that assembles the global nodal
% force vector
function [fitness] = CST_UDL_Fitness(q,t,weight,Ly,coords,elem,dofs,dofs_free,forces)

elements = size(elem,1);
fuerzas = zeros(length(forces),1);
count = 1;
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = elem(EL,1); n2 = elem(EL,2); % identify element node numbers
    n3 = elem(EL,3);
    x1 = coords(n1,1); y1 = coords(n1,2); % element node 1 - x,y coordinates
    x2 = coords(n2,1); y2 = coords(n2,2); % element node 2 - x,y coordinates
    x3 = coords(n3,1); y3 = coords(n3,2); % element node 3 - x,y coordinates
    
    dof12 = dofs(n1,2); % element node 1 - dofs
    dof22 = dofs(n2,2); % element node 2 - dofs
    dof32 = dofs(n3,2); % element node 3 - dofs
    
    % Nodal force vector
    x21 = x2 - x1; x31 = x3 - x1; % Triangle sides
    y21 = y2 - y1; y31 = y3 - y1;
    A = abs(x21*y31 - x31*y21)/2; % Area of element
    
    fuerzas(dof12) = fuerzas(dof12) + A*weight*t*1e-9/3; % Weight contribution
    fuerzas(dof22) = fuerzas(dof22) + A*weight*t*1e-9/3;
    fuerzas(dof32) = fuerzas(dof32) + A*weight*t*1e-9/3;
    
    if y1 == Ly  && y3 == Ly % UDL contribution
        fuerzas(dof12) = fuerzas(dof12) + q(count)*abs(x31)*1e-3/2;
        fuerzas(dof32) = fuerzas(dof32) + q(count)*abs(x31)*1e-3/2;
        count = count + 1;
        
    end
    
end

% Calculate norm ignoring restrained dofs (i.e. reactions)
fitness = forces(dofs_free) - fuerzas(dofs_free);


end