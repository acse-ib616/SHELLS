function [x_sol,fval,rel_tol,it] = GN_Truss(q,Ly,weight,A,elem,coords,dofs,dofs_free,forces,tol,maxIt,lambda)

nodes = length(coords)/2;
x_old = q;
rel_tol = 1;
it = 0;
% Calculate Jacobian
[j_values,j_rows,j_cols] = Truss_Jacobian(Ly,q,coords,elem);
j_rows = double(j_rows) + ones(size(j_rows));
j_cols = double(j_cols) + ones(size(j_cols));
jacobi = sparse(j_rows,j_cols,j_values,2*nodes,length(q));
jacob = jacobi(dofs_free,:);

while rel_tol > tol && it < maxIt
    
    [fval] = Truss_Fitness(x_old,weight,A,Ly,elem,coords,dofs,dofs_free,forces);
        
    p = (jacob'*jacob + lambda.*eye(length(x_old)))\(jacob'*fval);
    
    x_new = x_old + p;

    rel_tol = norm(fval);
    x_old = x_new;
    it = it + 1;
    

end

x_sol = x_new;
fval = norm(fval);

end