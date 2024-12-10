%% Draw Robot.
function plot_robot(handles, theta1, theta2, theta3, plotWorkspace, plotCoordinate)
    a2 = 0.65; 
    a3 = 0.65; 
    d1 = 0.25; 
    T1 = DHMatrix(theta1, d1, 0, pi/2);
    T2 = T1 * DHMatrix(theta2, 0, a2, 0);
    T3 = T2 * DHMatrix(theta3, 0, a3, 0);
    P0 = [0; 0; 0; 1];  
    P1 = T1 * P0;       
    P2 = T2 * P0;        
    P3 = T3 * P0;        
    axes(handles.axes1);
    cla;
    hold on;
    grid on;
    
    %% Draw link and base
    drawCylinderBetweenPoints(handles.axes1, P0, P1, 0.1, 0.1, [190/255, 40/255, 50/255]);
    drawCylinderBetweenPoints(handles.axes1, P1, P2, 0.07, 0.1, [22/255, 63/255, 155/255]);
    drawCylinderBetweenPoints(handles.axes1, P2, P3, 0.07, 0.1, [22/255, 63/255, 155/255]);
    %% Draw Joints
    drawHorizontalCylinder(handles.axes1, theta1, P1(1), P1(2), P1(3), 0.1, 0.3, [190/255, 40/255, 50/255]) ;
    drawHorizontalCylinder(handles.axes1, theta1, P2(1), P2(2), P2(3), 0.1, 0.3, [190/255, 40/255, 50/255]) ;

    %     %% Draw Tool
%     A = calculatePointOnLine(P2, P3, -0.25);
%     B = calculatePointOnLine(P2, P3, -0.251);
%     C = calculatePointOnLine(P2, P3, -0.1);
%     drawCylinderBetweenPoints(handles.axes1, A, B, 0.09, 0.5, [190/255, 40/255, 50/255]);
%     drawCylinderBetweenPoints(handles.axes1, B, P3, 0.02, 0.1, [190/255, 40/255, 50/255]);
%     drawCylinderBetweenPoints(handles.axes1, C, P3, 0.08, 0.1, [190/255, 40/255, 50/255]);
    
%     %% Draw Line on Cylinder
%     dir_P1P2 = P2 - P1;
%     length_P1P2 = norm(dir_P1P2);
%     t = linspace(-1, 2, 100); 
%     x = P1(1) + t * dir_P1P2(1);
%     y = P1(2) + t * dir_P1P2(2);
%     z = P1(3) + t * dir_P1P2(3);
%     plot3(x, y, z, 'LineWidth', 1.5, 'Color', 'black','LineStyle', '-  -');
%     %% draw Line on Cylinder
%     dir_P2P3 = P3 - P2;
%     length_P2P3 = norm(dir_P2P3);
%     x1 = P2(1) + t * dir_P2P3(1);
%     y1 = P2(2) + t * dir_P2P3(2);
%     z1 = P2(3) + t * dir_P2P3(3);
%     plot3(x1, y1, z1, 'LineWidth', 1.5, 'Color', 'black','LineStyle', '-  -');

    
    

    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis([-2 2 -2 2 -1 2]);  
    view(3);  
    rotate3d(handles.axes1, 'on'); 
    
    
    
  
 
    %% Draw coordinates 
    
    if plotCoordinate
    % At the base P0 (0, 0, 0)
    R0 = T1(1:3, 1:3);  % Rotation matrix for base (using T1)
    O0 = P0(1:3);  % Position at (0, 0, 0)
    
    % Draw the X, Y, Z axes at the base using the rotation matrix R0
    quiver3(O0(1), O0(2), O0(3), R0(1,1), R0(2,1), R0(3,1), 0.5, 'g', 'LineWidth', 2);  % Y-axis
    quiver3(O0(1), O0(2), O0(3), R0(1,2), R0(2,2), R0(3,2), 0.5, 'b', 'LineWidth', 2);  % Z-axis
    quiver3(O0(1), O0(2), O0(3), R0(1,3), R0(2,3), R0(3,3), 0.5, 'r', 'LineWidth', 2);  % X-axis
     % At P1
    R1 = T1(1:3, 1:3);  % Rotation matrix for P1
    O1 = P1(1:3);  % Position at P1
    quiver3(O1(1), O1(2), O1(3), R1(1,1), R1(2,1), R1(3,1), 0.5, 'r', 'LineWidth', 2);  % x-axis
    quiver3(O1(1), O1(2), O1(3), R1(1,2), R1(2,2), R1(3,2), 0.5, 'g', 'LineWidth', 2);  % y-axis
    quiver3(O1(1), O1(2), O1(3), R1(1,3), R1(2,3), R1(3,3), 0.5, 'b', 'LineWidth', 2);  % z-axis

    % At P2
    R2 = T2(1:3, 1:3);  % Rotation matrix for P2
    O2 = P2(1:3);  % Position at P2
    quiver3(O2(1), O2(2), O2(3), R2(1,1), R2(2,1), R2(3,1), 0.5, 'r', 'LineWidth', 2);  % x-axis
    quiver3(O2(1), O2(2), O2(3), R2(1,2), R2(2,2), R2(3,2), 0.5, 'g', 'LineWidth', 2);  % y-axis
    quiver3(O2(1), O2(2), O2(3), R2(1,3), R2(2,3), R2(3,3), 0.5, 'b', 'LineWidth', 2);  % z-axis
    end 
    
    % At P3
    R = T3(1:3, 1:3);  % Rotation matrix of the end-effector
    O = P3(1:3);  % End-effector position
    % Draw the X, Y, Z axes at the end-effector
    quiver3(O(1), O(2), O(3), R(1,1), R(2,1), R(3,1), 0.5, 'r', 'LineWidth', 2);  % x-axis
    quiver3(O(1), O(2), O(3), R(1,2), R(2,2), R(3,2), 0.5, 'g', 'LineWidth', 2);  % y-axis
    quiver3(O(1), O(2), O(3), R(1,3), R(2,3), R(3,3), 0.5, 'b', 'LineWidth', 2);  % z-axis
    
        %% Workspace 
%      if plotWorkspace
%         theta1_range = linspace(-180, 180, 15);
%         theta2_range = linspace(0, -180, 15);
%         theta3_range = linspace(-180, 0, 15); 
%         x_values = [];
%         y_values = [];
%         z_values = [];
%         for t1 = theta1_range
%             for t2 = theta2_range
%                 for t3 = theta3_range
%                     [x_workspace, y_workspace, z_workspace] = forward_kinematics(deg2rad(t1), deg2rad(t2), deg2rad(t3));             x_values = [x_values, x_workspace];
%                     y_values = [y_values, y_workspace];
%                     z_values = [z_values, z_workspace];
%                 end
%             end
%         end
%         
%         scatter3(x_values, y_values, z_values, 10, 'g', 'filled');
%     end
%     
    if plotWorkspace
    % Area 1
    theta1_range1 = linspace(-180, 180, 30);  
    theta2_range1 = linspace(-90, 90, 30);
    theta3_range1 = linspace(-90, 90, 30);
    
    % Area 2
    theta1_range2 = linspace(90,270, 30);  
    theta2_range2 = linspace(-10,-10,30);
%     theta2_range2 = linspace(0, 10, 20);
    theta3_range2 = linspace(-130, 25, 30);
    
    %Area 3
    theta1_range3 = linspace(-90, 90, 30);
    theta2_range3 = linspace(-10, -10, 30);
    theta3_range3 = linspace(-130, 25, 30);

    x_values = [];
    y_values = [];
    z_values = [];
    
  
    for t1 = theta1_range1
        for t2 = theta2_range1
            for t3 = theta3_range1
                [x_workspace, y_workspace, z_workspace] = forward_kinematics(deg2rad(t1), deg2rad(t2), deg2rad(t3));
                
                
                if z_workspace > 0
                    x_values = [x_values, x_workspace];
                    y_values = [y_values, y_workspace];
                    z_values = [z_values, z_workspace];
                end
            end
        end
    end
    
   
    scatter3(x_values, y_values, z_values, 10, 'g', 'filled');
    hold on;  
    
    
    x_values = [];
    y_values = [];
    z_values = [];
    
    for t1 = theta1_range2
        for t2 = theta2_range2
            for t3 = theta3_range2
                [x_workspace, y_workspace, z_workspace] = forward_kinematics(deg2rad(t1), deg2rad(t2), deg2rad(t3));
                
               
                if z_workspace < 0
                    x_values = [x_values, x_workspace];
                    y_values = [y_values, y_workspace];
                    z_values = [z_values, z_workspace];
                end
            end
        end
    end
    
 
    scatter3(x_values, y_values, z_values, 10, 'g', 'filled');
    
   
    x_values = [];
    y_values = [];
    z_values = [];
    
    for t1 = theta1_range3
        for t2 = theta2_range3
            for t3 = theta3_range3
                [x_workspace, y_workspace, z_workspace] = forward_kinematics(deg2rad(t1), deg2rad(t2), deg2rad(t3));
                
                x_values = [x_values, x_workspace];
                y_values = [y_values, y_workspace];
                z_values = [z_values, z_workspace];
            end
        end
    end

    scatter3(x_values, y_values, z_values, 10, 'g', 'filled');
    
    hold off; 
end
