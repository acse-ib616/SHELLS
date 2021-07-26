% This function calls the MEX C++ function that assembles the global nodal
% force vector
function fitness = CST_UDL_Fitness(q,t,weight,Ly,coords,elem,dofs_free,forces)

fuerzas = F_CST_FULL(t,weight,Ly,q,coords,elem);

% Calculate norm ignoring restrained dofs (i.e. reactions)
fitness = norm(fuerzas(dofs_free) - forces(dofs_free));

end