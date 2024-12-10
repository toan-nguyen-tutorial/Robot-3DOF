%% --- Inverse kinematics function
function [theta1, theta2, theta3] = inverse_kinematics(x, y, z)
    a2 = 0.65; 
    a3 = 0.65;  
    d1 = 0.25;  
    theta1 = atan2(y,x);
    r = sqrt(x^2 + y^2);
    c3 = (r^2 + (z - d1)^2 - a2^2 - a3^2) / (2 * a2 * a3); 
    s3 = sqrt(1 - c3^2); 
    theta3 = atan2(s3, c3);  
    c2 = (r * (a2 + a3 * c3) + (z - d1) * a3 * s3) / (a2^2 + a3^2 + 2 * a2 * a3 * c3);
    s2 = (-r * a3 * s3 + (z - d1) * (a2 + a3 * c3)) / (a2^2 + a3^2 + 2 * a2 * a3 * c3);
    theta2 = atan2(s2, c2);
