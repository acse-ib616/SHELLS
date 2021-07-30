function [jacob] = CST_UDL_jacobian(q,Ly,coords,elem,dofs,dofs_free)
elements = size(elem,1);
jacobi = zeros(length(dofs),length(q));
count = 1;
for EL = 1:elements % loop through all elements & build stiffness matrix
    n1 = elem(EL,1); % identify element node numbers
    n3 = elem(EL,3);
    x1 = coords(n1,1); y1 = coords(n1,2); % element node 1 - x,y coordinates
    x3 = coords(n3,1); y3 = coords(n3,2); % element node 3 - x,y coordinates
    
    dof12 = dofs(n1,2); % element node 1 - dofs
    dof32 = dofs(n3,2); % element node 3 - dofs
    x31 = x3 - x1;
    
    if y1 == Ly  && y3 == Ly % UDL contribution
        jacobi(dof12,count) = 0.5*abs(x31)*1e-3;
        jacobi(dof32,count) = 0.5*abs(x31)*1e-3;
        count = count + 1;
        
    end
    
end

jacob = jacobi(dofs_free,:);
end