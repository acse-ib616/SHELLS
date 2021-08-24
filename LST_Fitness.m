% This function calls the MEX C++ function that assembles the global nodal
% force vector
function [fitness] = LST_Fitness(q,t,weight,Ly,coords,elem,dofs,dofs_free,forces)

elements = length(elem)/6;
fuerzas = zeros(length(forces),1);
count = 1;
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = elem((EL-1)*6+1); n2 = elem((EL-1)*6+2); n3 = elem((EL-1)*6+3); % identify element node numbers
    n4 = elem((EL-1)*6+4); n5 = elem((EL-1)*6+5); n6 = elem((EL-1)*6+6);
    x1 = coords((n1-1)*2+1); y1 = coords((n1-1)*2+2); % element node 1 - x,y coordinates
    x2 = coords((n2-1)*2+1); y2 = coords((n2-1)*2+2); % element node 2 - x,y coordinates
    x3 = coords((n3-1)*2+1); y3 = coords((n3-1)*2+2); % element node 3 - x,y coordinates
    x4 = coords((n4-1)*2+1); y4 = coords((n4-1)*2+2); % element node 4 - x,y coordinates
    
    dof12 = dofs(n1,2); % element node 1 - dofs
    dof22 = dofs(n2,2); % element node 2 - dofs
    dof42 = dofs(n4,2); % element node 3 - dofs
    dof52 = dofs(n5,2); % element node 3 - dofs
    dof62 = dofs(n6,2); % element node 3 - dofs
    
    x21 = x2 - x1; x31 = x3 - x1; % Triangle sides
    y21 = y2 - y1; y31 = y3 - y1;
    A = abs(x21*y31 - x31*y21)/2; % Area of element
    
    x41 = x4 - x1;
    
    % Weight contribution. No vertex loads!
    fuerzas(dof42) = fuerzas(dof42) + A*weight*t*1e-9/3; 
    fuerzas(dof52) = fuerzas(dof52) + A*weight*t*1e-9/3;
    fuerzas(dof62) = fuerzas(dof62) + A*weight*t*1e-9/3;
    
    if y1 == Ly  && y2 == Ly && y4 == Ly % UDL contribution
        fuerzas(dof12) = fuerzas(dof12) + q(count)*abs(x41)*1e-3/2;
        fuerzas(dof22) = fuerzas(dof22) + q(count)*abs(x41)*1e-3/2;
        % Double contribution for midpoint along edge
        fuerzas(dof42) = fuerzas(dof42) + q(count)*abs(x41)*1e-3; 
        count = count + 1;
        
    end
    
end

% Calculate norm ignoring restrained dofs (i.e. reactions)
fitness = forces(dofs_free) - fuerzas(dofs_free);


end