% Iñigo Basterretxea Jacob 01246662

% This function assembles the finite element mesh for the truss, CST and
% LST elements

% Inputs:
% - Lx: width of domain
% - Ly: height of domain
% - nx: no. of divisions along x-direction
% - ny: no. of divisions along y-direction
% - element_type: CST, LST or truss element type
% - BC: boundary conditions

% Output:
% - NODES: struct with matrices of x,y nodal coordinates (NODES.coords) & DOFs (NODES.dofs)
% - ELEMENTS: matrix with element nodal connectivity
% - dofs_restrained: vector of restricted DOFs of system
% - dofs_free: vector of unrestricted DOFs of system
function [NODES,ELEMENTS,dofs_free,dofs_restrained] = mesh_FE(Lx,Ly,nx,ny,element_type,BC)

if element_type == 1 % CST
    elements = nx*2*ny; % no. of elements
    dx = Lx/nx; % Distance between nodes in x
    dy = Ly/ny; % Distance between nodes in y
    Nx = nx+1; % no. of nodes in x direction
    Ny = ny+1; % no. of nodes in y direction
    N = Nx*Ny; % no. of nodes

    % Note that node numbers are identified by their row number
    % Specifying nodal x-y coordinates. 1st column is x & 2nd is y
    NODES.coords = zeros(N,2);
    for i=1:Nx
        for j=1:Ny
            NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
        end
    end

    % Specifying  DOFs in pairs per node. Odd numbers are horizontal DOFs &
    % even numbers are vertical DOFs
    NODES.dofs = zeros(N,2);
    for i=1:N  
        NODES.dofs(i,:) = [2*i-1, 2*i];
    end

    % Specifying element nodal connectivity in rows of 3 nodes per element
    % (order of the elements does not matter)
    ELEMENTS = zeros(elements,3);
    for i=1:ny
        for j=1:nx
            ELEMENTS((i-1)*2*nx+2*j-1,:) = [i+(j-1)*(ny+1), i+ny+1+(j-1)*(ny+1), i+1+(j-1)*(ny+1)];
            ELEMENTS((i-1)*2*nx+2*j,:) = [i+(j-1)*(ny+1)+1, i+ny+1+(j-1)*(ny+1), i+ny+2+(j-1)*(ny+1)];
        end
    end
    
    if BC == 1
        % unknown nodal x-y dofs
        dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*Ny-2,2*Ny:2*N];
        % known nodal x-y dofs due to BC
        dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*Ny-1];
        
    else
        dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*N-Ny,2*N-Ny+2:2*N];
        dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*N-Ny+1];
    end

elseif element_type == 2 % LST
    
    elements = nx*2*ny; % no. of elements
    dx = Lx/(2*nx); % Distance between nodes in x
    dy = Ly/(2*ny); % Distance between nodes in y
    Nx = 2*nx+1; % no. of nodes in x direction
    Ny = 2*ny+1; % no. of nodes in y direction
    N = Nx*Ny; % no. of nodes
    
    % Note that node numbers are identified by their row number
    % Specifying nodal x-y coordinates. 1st column is x & 2nd is y
    NODES.coords = zeros(N,2);
    for i=1:Ny
        for j=1:Nx
            NODES.coords((i-1)*Nx+j,:) = [(j-1)*dx,dy*(i-1)];
        end
    end
    
    % Specifying  DOFs in pairs per node. Odd numbers are horizontal DOFs &
    % even numbers are vertical DOFs
    NODES.dofs = zeros(N,2);
    for i=1:N
        
        NODES.dofs(i,:) = [2*i-1, 2*i];
    end
    % Note that node numbers are identified by their row number
    
    % Specifying element nodal connectivity in rows of 3 nodes per element
    % (order of the elements does not matter)
    ELEMENTS = zeros(elements,6);
    for i=1:elements
        if mod(floor(i/nx),2) == 0 && mod(i,nx)~= 0
            ELEMENTS(i,:) = [2*mod(i,nx)-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0)+1, 2*Nx+2*mod(i,nx)-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), Nx+2*mod(i,nx)+Nx*floor(i/nx), Nx+2*mod(i,nx)-1+Nx*floor(i/nx)];
        elseif mod(i,nx)== 0 && mod(floor(i/nx),2) ~= 0
            ELEMENTS(i,:) = [Nx-2+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), Nx+2*Nx-2+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), Nx-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0), -1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0), -2+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),0)]; %Nx+2*Nx-2+2*Nx*max(floor(i/nx)-floor(i/(2*nx))-1,0)
        elseif mod(i,nx)== 0 && mod(floor(i/nx),2) == 0
            ELEMENTS(i,:) = [2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)+Nx, Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-2, Nx+2*Nx*(max(floor(i/nx)-floor(i/(2*nx)),1)-1), Nx-1+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1), Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx-1, Nx+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx];
        else
            ELEMENTS(i,:) = [2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)+2*mod(i,nx)+1, 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-1, 2*mod(i,nx)+1+2*Nx*(max(floor(i/nx)-floor(i/(2*nx)),1)-1), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1), 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx, 2*mod(i,nx)+2*Nx*max(floor(i/nx)-floor(i/(2*nx)),1)-Nx+1];
        end
    end
    
    if BC == 1
        % Degrees of freedom
        % unknown nodal x-y dofs
        dofs_free = [1:N-Nx,N-Nx+3:2*N-2*Nx,2*N-2*Nx+2:2*N];
        % known nodal x-y dofs due to BC
        dofs_restrained = [N-Nx+1,N-Nx+2,2*N-2*Nx+1];
        
    else
        dofs_free = [1:N-Nx,N-Nx+3:N+Nx-1,N+Nx+1:2*N];
        dofs_restrained = [N-Nx+1,N-Nx+2,N+Nx];
    end
    
else % Truss by default
    elements = nx*3*ny + nx + ny; % no. of elements
    dx = Lx/nx; % Distance between nodes in x
    dy = Ly/ny; % Distance between nodes in y
    Nx = nx+1; % Nodes in x direction
    Ny = ny+1; % Nodes in y direction
    N = Nx*Ny; % no. of nodes

    % Note that node numbers are identified by their row number
    % Specifying nodal x-y coordinates. 1st column is x & 2nd is y
    NODES.coords = zeros(N,2);
    for i=1:Nx
        for j=1:Ny
            NODES.coords((i-1)*Ny+j,:) = [(i-1)*dx,dy*(j-1)];
        end
    end
    
    % Specifying  DOFs in pairs per node. Odd numbers are horizontal DOFs &
    % even numbers are vertical DOFs
    NODES.dofs = zeros(N,2);
    for i=1:N
        NODES.dofs(i,:) = [2*i-1, 2*i];
    end

    % Specifying element nodal connectivity in rows of 2 nodes per element
    % (order of the elements does not matter)
    ELEMENTS = zeros(elements,2);
    for i = 0:ny-1
        % Horizontal member
        for k = 1:nx
            ELEMENTS((3*nx+1)*i+k,:) = [i+1, i+ny+2] + [(k-1)*(ny+1), (k-1)*(ny+1)]; 
        end
        % Vertical member
        for j = 1:nx+1
            ELEMENTS((3*nx+1)*i+nx+j,:) = [i+1, i+2] + [(j-1)*(ny+1), (j-1)*(ny+1)]; 
        end
        % High-low diagonal
        for l = 1:nx
            ELEMENTS((3*nx+1)*i+2*nx+1+l,:) = [i+2, i+ny+2] + [(l-1)*(ny+1), (l-1)*(ny+1)]; 
        end
    end
    % Top horizontal members
    for k = 1:nx
            ELEMENTS((3*nx+1)*ny+k,:) = [ny+1, ny+ny+2] + [(k-1)*(ny+1), (k-1)*(ny+1)];
    end
    
    if BC == 1
        % Degrees of freedom & other parameters
        % unknown nodal x-y dofs
        dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*Ny-2,2*Ny:2*N]; 
        % known nodal x-y dofs due to BC
        dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*Ny-1];
    else
        dofs_free = [1:(floor(ny/2)+1)*2-2,(floor(ny/2)+1)*2+1:2*N-Ny,2*N-Ny+2:2*N];
        dofs_restrained = [(floor(ny/2)+1)*2-1,(floor(ny/2)+1)*2,2*N-Ny+1];
    end
  
end

end