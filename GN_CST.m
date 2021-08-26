function [x_sol,fval,it] = GN_CST(q,t,weight,Ly,coords,elem,dofs,dofs_free,forces,tol,maxIt,lambda)
% This function returns the value of the UDL applied at each of the top elements,
% using the Gauss-Newton method with the Levenberg–Marquardt modification


nodes = length(coords)/2; % no. of nodes

x_old = q; % Initialise previous solution estimate
rel_tol = 1; % Initialise tolerance
it = 0; % Initialise iteration counter

% Calculate sparse Jacobian for CST element by calling C++ function
[j_values,j_rows,j_cols] = CST_Jacobian(Ly,q,coords,elem);
% Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
% position and column index vectors
j_rows = double(j_rows) + ones(size(j_rows));
j_cols = double(j_cols) + ones(size(j_cols));
jacobi = sparse(j_rows,j_cols,j_values,2*nodes,length(q));

% Ignore restricted DOFs
jacob = jacobi(dofs_free,:);

while rel_tol > tol && it < maxIt
    
    % Evaluate function value (i.e. the residual R = y - f(x))
    fval = CST_UDL_Fitness(x_old,t,weight,Ly,coords,elem,dofs,dofs_free,forces);
    
    % Calculate step
    p = (jacob'*jacob + lambda.*eye(length(x_old)))\(jacob'*fval);
    
    % Update solution estimate
    x_new = x_old + p;

    % Calculate L2-norm
    rel_tol = norm(fval);
    
    % Update previous iteration estimate for next iteration
    x_old = x_new;
    
    % Update iteration counter
    it = it + 1;
    

end

% Output solution
x_sol = x_new;

% Output L2-norm
fval = norm(fval);

end