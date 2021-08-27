% Iñigo Basterretxea Jacob 01246662

% This function solves the 2 matrix equations to obtain the nodal
% displacements and the nodal forces vectors U & F

% Inputs:
% - K: sparse global stiffness matrix
% - dofs_restrained: vector of restricted DOFs of system
% - dofs_free: vector of unrestricted DOFs of system
% - forces: external nodal forces vector

% Output:
% - U: struct with matrices of x,y nodal coordinates (NODES.coords) & DOFs (NODES.dofs)
% - F: nodal forces vector including reactions
function [U,F] = solve_FE(K,forces,dofs_free,dofs_restrained)
F = forces;
N = length(forces)/2; % no. of nodes

% Specification of submatrices directly
KRR = K(dofs_restrained,dofs_restrained);
KRF = K(dofs_restrained,dofs_free);
KFR = K(dofs_free,dofs_restrained);
KFF = K(dofs_free,dofs_free);
fF = F(dofs_free);
% Cannot yet form fR as the reactions are unknown

% Solution for the unknown nodal dofs
uR = zeros(length(dofs_restrained),1); % BC - zero displacement on restricted nodes
uF = KFF\(fF - KFR*uR); % 1st matrix equation
U = zeros(2*N,1);
U(dofs_restrained) = uR;
U(dofs_free) = uF; % full nodal dof vector

% Solution for the unknown reactions
fR = KRF*uF + KRR*uR; % 2nd matrix equation
F(dofs_restrained) = fR; % full nodal force vector

% Test force balance
if abs(sum(F)/max(abs(F)))*100 > 1e-4
    msg = 'ABORTED: there is not force balance';
    error(msg);
end
end