function [fitness] = Truss_Fitness(q,weight,A,Ly,elem,coords,dofs,dofs_free,forces)

elements = length(elem)/2;
fuerzas = zeros(length(forces),1);
count = 1;
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = elem((EL-1)*2+1); n2 = elem((EL-1)*2+2); % identify element node numbers
    dof12 = dofs(n1,2); % element node 1 - dofs
    dof22 = dofs(n2,2); % element node 2 - dofs
    x1 = coords((n1-1)*2+1); y1 = coords((n1-1)*2+2); % element node 1 - x,y coordinates
    x2 = coords((n2-1)*2+1); y2 = coords((n2-1)*2+2); % element node 2 - x,y coordinates
    
    x21 = x2 - x1;
    
    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
    
    fuerzas(dof12) = fuerzas(dof12) + A*weight*Le*1e-9/3; % Weight contribution
    fuerzas(dof22) = fuerzas(dof22) + A*weight*Le*1e-9/3;
    
    if y1 == Ly  && y2 == Ly % UDL contribution
        fuerzas(dof12) = fuerzas(dof12) + q(count)*abs(x21)*1e-3/2;
        fuerzas(dof22) = fuerzas(dof22) + q(count)*abs(x21)*1e-3/2;
        count = count + 1;
    end
    
end


% Calculate norm ignoring restrained dofs (i.e. reactions)
fitness =  forces(dofs_free) - fuerzas(dofs_free);

end