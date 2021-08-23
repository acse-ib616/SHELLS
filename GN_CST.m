function [x_sol,fval,rel_tol,it] = GN_CST(q,t,weight,Ly,coords,elem,dofs,dofs_free,forces,tol,maxIt,lambda)
% This function returns the value of the UDL applied at each of the top elements,
% using the Gauss-Newton method with the Levenberg–Marquardt modification



x_old = q; % Solution at previous iteration
rel_tol = 1; 
it = 0; % Iteration counter

% Calculate Jacobian
jacobi = CST_Jacobian(Ly,q,coords,elem);
jacob = jacobi(dofs_free,:);

while rel_tol > tol && it < maxIt
    
    % Evaluate function value (i.e. the residual R = y - f(x))
    fval = CST_UDL_Fitness(x_old,t,weight,Ly,coords,elem,dofs,dofs_free,forces);
    
    p = (jacob'*jacob + lambda.*eye(length(x_old)))\(jacob'*fval);
    
    % Update solution
    x_new = x_old + p;

    rel_tol = norm(fval); % Euclidean norm of R = y - f(x)
    x_old = x_new; % Update previous iteration solution with new solution
    it = it + 1; % Increase iteration counter
    

end

x_sol = x_new;

% Fitness value is the Euclidean norm of R = y - f(x)
fval = norm(fval);

end