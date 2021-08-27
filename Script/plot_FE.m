% Iñigo Basterretxea Jacob 01246662
% This function is largely based on code written by Dr. Adam Sadowski
function p = plot_FE(ELEMENTS,COORDS,DOFS,U,amp,element_type,solve_U)

elements = size(ELEMENTS,1); % no. of nodes

if solve_U
    % Updating nodal coordinates after deformation
    % Creating vector with amplified deformed coordinate for plotting purposes
    amp_coords = zeros(size(COORDS));
    for I = 1:size(COORDS,1)
        for J = 1:size(COORDS,2)
            amp_coords(I,J) = COORDS(I,J) + U(DOFS(I,J))*amp;
        end
    end
end

% Plotting
p = figure('units','normalized','outerposition',[0 0 1 1]);
hold all; grid on;
axis equal;

if element_type == 1 || element_type == 2 % CST or LST
    for EL = 1:elements
        n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); n3 = ELEMENTS(EL,3); % identify element node numbers
        
        % Plotting original structure
        x1 = COORDS(n1,1); y1 = COORDS(n1,2); % element node 1 - x,y original coordinates
        x2 = COORDS(n2,1); y2 = COORDS(n2,2); % element node 2 - x,y original coordinates
        x3 = COORDS(n3,1); y3 = COORDS(n3,2); % element node 3 - x,y original coordinates
        patch([x1,x2,x3],[y1,y2,y3],[0.5 0.5 0.5]);
        
        if solve_U
            % Plotting amplified deformed structure
            x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
            x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
            x3_amp = amp_coords(n3,1); y3_amp = amp_coords(n3,2); % element node 3 - x,y amplified deformed coordinates
            patch([x1_amp,x2_amp,x3_amp],[y1_amp,y2_amp,y3_amp],'r');
        end
    end
    
else % Truss by default
    for EL = 1:elements
        % identify element node numbers
        n1 = ELEMENTS(EL,1); n2 = ELEMENTS(EL,2); 
        
        % Plotting original structure
        x1 = COORDS(n1,1); y1 = COORDS(n1,2); % element node 1 - x,y original coordinates
        x2 = COORDS(n2,1); y2 = COORDS(n2,2); % element node 2 - x,y original coordinates
        plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3);
        
        if solve_U
            % Plotting amplified deformed structure
            x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
            x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3);
        end
    
    end
end
xlabel('x coordinate');
ylabel('y coordinate');
if solve_U
    title('Deformed shape');
else
    title('Mesh');
end
set(gca,'FontSize',30);
end