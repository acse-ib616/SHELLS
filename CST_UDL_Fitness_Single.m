% This function calls the MEX C++ function that assembles the global nodal
% force vector
function fitness = CST_UDL_Fitness_Single(q,t,weight,Ly,coords,elem,dofs_free,forces)

fuerzas = F_CST_UDL(t,q,weight,Ly,coords,elem);

% Calculate norm ignoring restrained dofs (i.e. reactions)
fitness = norm(fuerzas(dofs_free) - forces(dofs_free));

end