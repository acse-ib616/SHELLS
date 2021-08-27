% Iñigo Basterretxea Jacob 01246662

% This function assembles the global stiffness matrix K sparsely by calling
% a C++ MEX function for the truss, CST and LST elements

% Inputs:
% - COORDS: vector of x,y nodal coordinates
% - ELEMENTS: vector of element nodal connectivity
% - E: Young's modulus
% - t: element thickness or cross-section if it is a truss element
% - nu: Poisson ratio
% - element_type: CST, LST or truss element type

% Output:
% - K: sparse global stiffness matrix
function [K] = K_FE(ELEMENTS,COORDS,E,t,nu,element_type)

% Make coords and elements matrices 1-dimensional. Numbers are double by
% default in MATLAB, so we must convert to unsigned integer
coords = reshape(COORDS',1,[]);
elem = uint32(reshape(ELEMENTS',1,[]));

N = size(COORDS,1); % no. of nodes

% Call C++ global stiffness matrix assembler
[k_values,k_rows,k_cols]  = FE_K(E,nu,t,coords,elem,uint32(element_type));
% Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
% position and column index vectors
k_rows = double(k_rows) + ones(size(k_rows));
k_cols = double(k_cols) + ones(size(k_cols));
K = sparse(k_rows,k_cols,k_values,2*N,2*N);


% Test K symmetry
if max(max(K-K.')) > 1e-8
    msg = 'ABORTED: stiffness matrix K is NOT symmetric.';
    error(msg);
end
end