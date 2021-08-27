% Iñigo Basterretxea Jacob 01246662
function [K] = K_FE(ELEMENTS,COORDS,E,t,nu,element_type)

% Make coords and elements matrices 1-dimensional. Numbers are double by
% default in MATLAB, so we must convert to unsigned integer
coords = reshape(COORDS',1,[]);
elem = uint32(reshape(ELEMENTS',1,[]));

N = size(COORDS,1); % no. of nodes

if element_type == 1 % CST
    % Call C++ LST stiffness matrix assembler
    [k_values,k_rows,k_cols]  = CST_K(E,nu,t,coords,elem);
    % Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
    % position and column index vectors
    k_rows = double(k_rows) + ones(size(k_rows));
    k_cols = double(k_cols) + ones(size(k_cols));
    K = sparse(k_rows,k_cols,k_values,2*N,2*N);
    

elseif element_type == 2 % LST
    % Call C++ LST stiffness matrix assembler
    [k_values,k_rows,k_cols]  = LST_K(E,nu,t,coords,elem);
    % Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
    % position and column index vectors
    k_rows = double(k_rows) + ones(size(k_rows));
    k_cols = double(k_cols) + ones(size(k_cols));
    K = sparse(k_rows,k_cols,k_values,2*N,2*N);

else % Truss by default
    [k_values,k_rows,k_cols] = Truss_K(E,t,coords,elem);
    % Matlab is 1-index based. C++ is 0-index based, so we must add 1 to row
    % position and column index vectors
    k_rows = double(k_rows) + ones(size(k_rows));
    k_cols = double(k_cols) + ones(size(k_cols));
    K = sparse(k_rows,k_cols,k_values,2*N,2*N);
end

end