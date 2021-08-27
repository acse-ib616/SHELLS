% Iñigo Basterretxea Jacob 01246662

% This function performs the Gauss-Newton method to optimise for a vector
% of top boundary UDLs a rectangular system of finite elements

% Inputs:
% - x0: initial estimate for vector of top boundary UDLs
% - options: Gauss-Newton optimisation parametres (eg. tol,maxIt,lambda)
% - COORDS: vector of x,y nodal coordinates
% - ELEMENTS: vector of element nodal connectivity
% - DOFS: vector of DOFs of system
% - dofs_free: vector of unrestricted DOFs of system
% - forces: input nodal forces vector F
% - constants: element constants (eg. E,t,weight)
% - element_type: CST, LST or truss element type

% Outputs:
% - xsol: vector of solutions for q
% - fval: L2-norm of input F vs assembled F
% - it: no. of iterations performed
function [x_sol,fval,it] = load_reconstruction(x0,forces,ELEMENTS,COORDS,DOFS,dofs_free,options,constants,element_type)
tol = options(1); % input L2-norm tolerance
maxIt = options(2); % maximum no. of iterations
lambda = options(3); % Levenberg–Marquardt Modification coefficient
Ly = constants(1); % height of domain
t = constants(2); % thickness of element (CST,LST) or cross-sectional area (truss)
weight = constants(3); % unit weight of material

% Make coords and elements matrices 1-dimensional. Numbers are double by
% default in MATLAB, so we must convert to unsigned integer
coords = reshape(COORDS',1,[]);
elem = uint32(reshape(ELEMENTS',1,[]));

N = size(COORDS,1); % no. of nodes
x_old = x0; % Initialise previous solution estimate
rel_tol = 1; % Initialise tolerance
it = 0; % Initialise iteration counter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% ASSEMBLE JACOBIAN %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate sparse Jacobian by calling C++ function
[j_values,j_rows,j_cols] = FE_Jacobian(Ly,x_old,coords,elem,uint32(element_type));
% Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
% position and column index vectors
j_rows = double(j_rows) + ones(size(j_rows));
j_cols = double(j_cols) + ones(size(j_cols));
jacobi = sparse(j_rows,j_cols,j_values,2*N,length(x_old));
% Ignore restricted DOFs
jacob = jacobi(dofs_free,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% GAUSS-NEWTON %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while rel_tol > tol && it < maxIt
    
    % Evaluate function value (i.e. the residual R = y - f(x))
    F = F_FE(ELEMENTS,COORDS,DOFS,Ly,t,x_old,weight,element_type);
    % Calculate residual ignoring restricted DOFs (i.e. reactions)
    fval = forces(dofs_free) - F(dofs_free);
    
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

if it == maxIt && rel_tol > tol
    fprintf('WARNING: did not after maximum no. of iterations \n');
end

% Output solution
x_sol = x_new;

% Output L2-norm
fval = norm(fval);
end