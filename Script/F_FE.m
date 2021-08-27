% Iñigo Basterretxea Jacob 01246662
function [F] = F_FE(ELEMENTS,COORDS,DOFS,Ly,t,q,weight,element_type)

N = size(COORDS,1); % no. of nodes
elements = size(ELEMENTS,1); % no. of elements
counter = 1; % Initialise q-vector counter
F = zeros(2*N,1); % Initialise nodal forces vector

if element_type == 1 % CST
    % Constructing the nodal forces vector
    for EL = 1:elements % loop through all elements & build stiffness matrix
        % identify element node numbers
        n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3);
        
        % identify node coordinates
        x1 = COORDS(n1,1); y1 = COORDS(n1,2);
        x2 = COORDS(n2,1); y2 = COORDS(n2,2);
        x3 = COORDS(n3,1); y3 = COORDS(n3,2);
        
        % identify element DOFs
        dof12 = DOFS(n1,2);
        dof22 = DOFS(n2,2);
        dof32 = DOFS(n3,2);
        
        % Triangle sides
        x21 = x2 - x1; x31 = x3 - x1; y21 = y2 - y1; y31 = y3 - y1;
        
        % Area of element
        A = abs(x21*y31 - x31*y21)/2;
        
        % Weight contribution at all elements
        F(dof12) = F(dof12) + A*weight*t*1e-9/3;
        F(dof22) = F(dof22) + A*weight*t*1e-9/3;
        F(dof32) = F(dof32) + A*weight*t*1e-9/3;
        
        % UDL contribution at top boundary elements
        if y1 == Ly  && y3 == Ly
            F(dof12) = F(dof12) + q(counter)*abs(x31)*1e-3/2;
            F(dof32) = F(dof32) + q(counter)*abs(x31)*1e-3/2;
            counter = counter + 1;
        end
        
    end
    

elseif element_type == 2 % LST
    % Constructing the nodal forces vector
    for EL = 1:elements % loop through all elements & build stiffness matrix
        
        % identify element node numbers
        n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3);
        n4 = ELEMENTS(EL,4); n5 = ELEMENTS(EL,5); n6 = ELEMENTS(EL,6);
        
        % identify node coordinates
        x1 = COORDS(n1,1); y1 = COORDS(n1,2);
        x2 = COORDS(n2,1); y2 = COORDS(n2,2);
        x3 = COORDS(n3,1); y3 = COORDS(n3,2);
        x4 = COORDS(n4,1); y4 = COORDS(n4,2);
        
        % identify element DOFs
        dof12 = DOFS(n1,2);
        dof22 = DOFS(n2,2);
        dof42 = DOFS(n4,2);
        dof52 = DOFS(n5,2);
        dof62 = DOFS(n6,2);
        
        % Triangle sides
        x21 = x2 - x1; x31 = x3 - x1;
        y21 = y2 - y1; y31 = y3 - y1;
        
        % Half of a side
        x41 = x4 - x1;
        
        % Area of element
        A = abs(x21*y31 - x31*y21)/2;
        
        % Weight contribution at all elements. No vertex loads!
        F(dof42) = F(dof42) + A*weight*t*1e-9/3;
        F(dof52) = F(dof52) + A*weight*t*1e-9/3;
        F(dof62) = F(dof62) + A*weight*t*1e-9/3;
        
        % UDL contribution at top boundary elements
        if y1 == Ly  && y2 == Ly && y4 == Ly
            F(dof12) = F(dof12) + q(counter)*abs(x41)*1e-3/2;
            F(dof22) = F(dof22) + q(counter)*abs(x41)*1e-3/2;
            
            % Double contribution for midpoint along edge
            F(dof42) = F(dof42) + q(counter)*abs(x41)*1e-3;
            counter = counter + 1;
            
        end
        
    end


else % Truss by default
    % Constructing the nodal forces vector
    for EL = 1:elements % loop through all elements & build stiffness matrix
        
        % identify element node numbers
        n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2);
        
        % identify node coordinates
        x1 = COORDS(n1,1); y1 = COORDS(n1,2);
        x2 = COORDS(n2,1); y2 = COORDS(n2,2);
        
        % identify element DOFs
        dof12 = DOFS(n1,2);
        dof22 = DOFS(n2,2);
        
        % Element length
        Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 );
        
        % Length projected in x
        x21 = x2 - x1;
        
        % Weight contribution at all elements
        F(dof12) = F(dof12) + t*weight*Le*1e-9/3;
        F(dof22) = F(dof22) + t*weight*Le*1e-9/3;
        
        % UDL contribution at top boundary elements
        if y1 == Ly  && y2 == Ly
            F(dof12) = F(dof12) + q(counter)*abs(x21)*1e-3/2;
            F(dof22) = F(dof22) + q(counter)*abs(x21)*1e-3/2;
            counter = counter + 1;
        end
        
    end
    
end
end