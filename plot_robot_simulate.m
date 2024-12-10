%% Draw Robot.

function plot_robot_simulate(handles, theta1, theta2, theta3)
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
    drawHorizontalCylinder(handles.axes1,theta1, P1(1), P1(2), P1(3), 0.1, 0.3, [190/255, 40/255, 50/255]) 
    drawHorizontalCylinder(handles.axes1,theta1, P2(1), P2(2), P2(3), 0.1, 0.3, [190/255, 40/255, 50/255])
    %% Draw Tool
    A = calculatePointOnLine(P2, P3, -0.25);
    B = calculatePointOnLine(P2, P3, -0.251);
    C = calculatePointOnLine(P2, P3, -0.1);
%     drawCylinderBetweenPoints(handles.axes1, A, B, 0.09, 0.5, [190/255, 40/255, 50/255]);
%     drawCylinderBetweenPoints(handles.axes1, B, P3, 0.02, 0.1, [190/255, 40/255, 50/255]);
%     drawCylinderBetweenPoints(handles.axes1, C, P3, 0.08, 0.1, [190/255, 40/255, 50/255]);
    
    %% Draw coordinates
    
    % At P3
    R = T3(1:3, 1:3);  % Rotation matrix of the end-effector
    O = P3(1:3);  % End-effector position
    % Draw the X, Y, Z axes at the end-effector
    quiver3(O(1), O(2), O(3), R(1,1), R(2,1), R(3,1), 0.5, 'r', 'LineWidth', 2);  % Y-axis
    quiver3(O(1), O(2), O(3), R(1,2), R(2,2), R(3,2), 0.5, 'g', 'LineWidth', 2);  % Z-axis
    quiver3(O(1), O(2), O(3), R(1,3), R(2,3), R(3,3), 0.5, 'b', 'LineWidth', 2);  % X-axis
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

    
  
    