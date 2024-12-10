%% --- Forward kinematics function
function [x, y, z] = forward_kinematics(theta1, theta2, theta3)
    a2 = 0.65; 
    a3 = 0.65; 
    d1 = 0.25;  % Lengths of links and offsets
    % Compute transformation matrices for each joint using DH convention
    T1 = DHMatrix(theta1, d1, 0, pi/2);  % Transformation from base to joint 1
    T2 = T1 * DHMatrix(theta2, 0, a2, 0);  % Transformation from joint 1 to joint 2
    T3 = T2 * DHMatrix(theta3, 0, a3, 0);  % Transformation from joint 2 to joint 3
    % Final velocity1 of the end-effector (P) in homogeneous coordinates
    P = T3 * [0; 0; 0; 1];  % Multiply with [0; 0; 0; 1] to get position
    % Extract the x, y, z components from the resulting vector
    x = P(1); 
    y = P(2); 
    z = P(3);
