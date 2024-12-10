function varargout = GUIDE(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function GUIDE_OpeningFcn(hObject, eventdata, handles, varargin)
axes(handles.axes1);
cla;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([-2 2 -2 2 -0.5 2]);
view(3);
%% Open Default (theta1;theta2;theta3)
plot_robot_simulate(handles, 0, 0, 0);
handles.output = hObject;
guidata(hObject, handles);
function varargout = GUIDE_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
function theta1_Callback(hObject, eventdata, handles)
function theta1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function theta2_Callback(hObject, eventdata, handles)
function theta2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function theta3_Callback(hObject, eventdata, handles)
function theta3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_x_Callback(hObject, eventdata, handles)
function edit_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_y_Callback(hObject, eventdata, handles)
function edit_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor'   ,'white');
end
function edit_z_Callback(hObject, eventdata, handles)
function edit_z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Button Foward Kinematics
function button_fwd_Callback(hObject, eventdata, handles)
theta1 = str2double(get(handles.theta1, 'String'));
theta2 = str2double(get(handles.theta2, 'String'));
theta3 = str2double(get(handles.theta3, 'String'));
if isnan(theta1) || isnan(theta2) || isnan(theta3)
   msgbox('Please enter valid angles for Theta1, Theta2, and Theta3.', 'Error', 'error');
    return;
end
% Convert Degree to Radian
theta1_rad = deg2rad(theta1); 
theta2_rad = deg2rad(theta2);
theta3_rad = deg2rad(theta3);    
%Draw Robot
plotCoordinate = handles.plotCoordinate;
plotWorkspace = handles.plotWorkspace; %call check box
plot_robot(handles, theta1_rad, theta2_rad, theta3_rad, plotWorkspace, plotCoordinate); 
% Forward kinematics to calculate X, Y, Z with theta is radians
[x, y, z] = forward_kinematics(theta1_rad, theta2_rad, theta3_rad);
% Display XYZ in GUI
set(handles.edit_x, 'String', sprintf('%.6f', x));
set(handles.edit_y, 'String', sprintf('%.6f', y));
set(handles.edit_z, 'String', sprintf('%.6f', z));
  
%% Button Inverse Kinematics
function button_inv_Callback(hObject, eventdata, handles)
x = str2double(get(handles.edit_x, 'String'));
y = str2double(get(handles.edit_y, 'String'));
z = str2double(get(handles.edit_z, 'String'));
if isnan(x) || isnan(y) || isnan(z)
    msgbox('Please enter valid coordinates for X, Y, and Z.', 'Error', 'error');
    return;
end
try
 % Inverse kinematics calculate theta1, theta2, theta3 
 [theta1, theta2, theta3] = inverse_kinematics(x, y, z);
    theta1_deg = rad2deg(theta1);
    theta2_deg = rad2deg(theta2);
    theta3_deg = rad2deg(theta3);
 % Display theta1, theta2, theta3 in GUI
   set(handles.theta1, 'String', sprintf('%.6f', theta1_deg));
   set(handles.theta2, 'String', sprintf('%.6f', theta2_deg));
   set(handles.theta3, 'String', sprintf('%.6f', theta3_deg));   
 % Draw Robot
   plotCoordinate = handles.plotCoordinate;
   plotWorkspace = handles.plotWorkspace; % call check box
   plot_robot(handles, theta1, theta2, theta3, plotWorkspace, plotCoordinate);
catch
   msgbox('Position is not reachable.', 'Error', 'error');
end


%% BUTTON HOME use Forward theta1 = 0; theta2 = 120; theta3 = -120
function btn_Home_Callback(hObject, eventdata, handles)
set(handles.theta1, 'String', '0');
set(handles.theta2, 'String', '120');
set(handles.theta3, 'String', '-120');
set(handles.edit_x, 'String', '0.00');
set(handles.edit_y, 'String', '0.00');
set(handles.edit_z, 'String', '0.00');
theta1 = str2double(get(handles.theta1, 'String'));
theta2 = str2double(get(handles.theta2, 'String'));
theta3 = str2double(get(handles.theta3, 'String'));
if isnan(theta1) || isnan(theta2) || isnan(theta3)
    msgbox('Please enter valid angles for Theta1, Theta2, and Theta3.', 'Error', 'error');
    return;
end


%Convert degree to radian
theta1_rad = deg2rad(theta1);
theta2_rad = deg2rad(theta2);
theta3_rad = deg2rad(theta3);

plotCoordinate = handles.plotCoordinate;
% Draw robot with theta is radian
plotWorkspace = handles.plotWorkspace; 
% G?i hàm v? robot v?i các giá tr? góc ?ã nh?p và tr?ng thái checkbox
plot_robot(handles, theta1_rad, theta2_rad, theta3_rad, plotWorkspace, plotCoordinate); 

%% Motion from A to B 
function btnMove_Callback(hObject, eventdata, handles)
    thetaA_deg = [45, 30, -30];  % Start point A (theta1, theta2, theta3) in degrees
    thetaB_deg = [-15, 30, 30];  % End point B (theta1, theta2, theta3) in degrees
    thetaM_deg = [34.33, -4.68, 57.68];
    thetaN_deg = [2.8, 12.136, 57.74];
    %Convert Degree to Radian
    thetaA = deg2rad(thetaA_deg);
    thetaB = deg2rad(thetaB_deg);
    thetaM = deg2rad(thetaM_deg);
    thetaN = deg2rad(thetaN_deg);
    %Calculate xyz use forward kinematics
    [x1, y1, z1] = forward_kinematics(thetaA(1), thetaA(2), thetaA(3));
    [x2, y2, z2] = forward_kinematics(thetaB(1), thetaB(2), thetaB(3));
    [x3, y3, z3] = forward_kinematics(thetaM(1), thetaM(2), thetaM(3));
    [x4, y4, z4] = forward_kinematics(thetaN(1), thetaN(2), thetaN(3));
    
    % Calculate distance
    d = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2); 
    % Parameter Vehicity max and ACC max
    v_max = 0.175; 
    a_max = 0.05;  
    % Calculate Time ACC and DEC
    t_acc = sqrt(d / (2 * a_max));  %use formulas "d = 0.5 * a * t^2"
    t_dec = t_acc; 
    % Calculate distance AM:=> ACC ;  MN:=> CONST ; NB:=>DEC
    d_acc = 0.5 * a_max * t_acc^2; % d_AM
    d_dec = 0.5 * a_max * t_dec^2; % d_NB
    d_constant = d - d_acc - d_dec; % d_MN
    if d_constant > 0
        t_constant = d_constant / v_max; 
    else
        t_constant = 0;
    end
    % Total time from A to B
    t_total = t_acc + t_constant + t_dec;
    % step for draw
    num_steps = 100; 
   
    time_acc = linspace(0, t_acc, num_steps * (t_acc / t_total));  
    time_const = linspace(t_acc, t_acc + t_constant, num_steps * (t_constant / t_total)); 
    time_dec = linspace(t_acc + t_constant, t_total, num_steps * (t_dec / t_total));  
    
    
    v_acc = a_max * time_acc; 
    v_const = v_max * ones(1, length(time_const));  
    v_dec = v_max - a_max * (time_dec - (t_acc + t_constant));  
    %% Draw Velocity of Robot
    time = [time_acc, time_const, time_dec];  
    velocity = [v_acc, v_const, v_dec];  
    axes(handles.velocity);  
    plot(time, velocity, 'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('Velocity');
    title('Velocity');
    grid on;
    
    
    %% Draw Position 
    pos_acc = 0.5 * a_max * time_acc.^2; 
    pos_const = d_acc + v_max * (time_const - t_acc);  
    pos_dec = d_acc + d_constant + v_max * (time_dec - (t_acc + t_constant)) - 0.5 * a_max * (time_dec - (t_acc + t_constant)).^2; 
    position = [pos_acc, pos_const, pos_dec];
    axes(handles.position);  
    plot(time, position, 'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('Position');
    title('Position');
    grid on;
    
    %% Draw Acceleration
    acceleration_acc = a_max * ones(1, length(time_acc));
    acceleration_const = zeros(1, length(time_const));
    acceleration_dec =  -a_max * ones(1, length(time_dec));
    acceleration = [acceleration_acc, acceleration_const, acceleration_dec];
    
    axes(handles.acceleration);  
    plot(time, acceleration, 'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('acceleration');
    title('Acceleration');
    grid on;
    
    %% Time steps
    num_steps_acc = round(num_steps * (t_acc / t_total));        
    num_steps_const = round(num_steps * (t_constant / t_total));  
    num_steps_dec = num_steps - num_steps_acc - num_steps_const; 

    %% ACC 
    theta1_acc = linspace(thetaA(1), thetaM(1), num_steps_acc);
    theta2_acc = linspace(thetaA(2), thetaM(2), num_steps_acc);
    theta3_acc = linspace(thetaA(3), thetaM(3), num_steps_acc);
    %% CONST
    theta1_const = linspace(thetaM(1), thetaN(1), num_steps_const);  
    theta2_const = linspace(thetaM(2), thetaN(2), num_steps_const);
    theta3_const = linspace(thetaM(3), thetaN(3), num_steps_const);
    %% DEC
    theta1_dec = linspace(thetaN(1), thetaB(1), num_steps_dec); 
    theta2_dec = linspace(thetaN(2), thetaB(2), num_steps_dec);
    theta3_dec = linspace(thetaN(3), thetaB(3), num_steps_dec);

    %% Final
    theta1 = [theta1_acc, theta1_const, theta1_dec];
    theta2 = [theta2_acc, theta2_const, theta2_dec];
    theta3 = [theta3_acc, theta3_const, theta3_dec];
    
    trajectory_x = zeros(1, length(time));
    trajectory_y = zeros(1, length(time));
    trajectory_z = zeros(1, length(time));
   %% Kh?i t?o m?ng
    theta1_dot_array = zeros(1, length(time));
    theta2_dot_array = zeros(1, length(time));
    theta3_dot_array = zeros(1, length(time));
        
    %% V? Theta
    axes(handles.plotTheta1);
    theta1_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta1 (degrees)');
    title('Theta1');
    grid on;

    axes(handles.plotTheta2);
    theta2_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta2 (degrees)');
    title('Theta2');
    grid on;
  
    axes(handles.plotTheta3);
    theta3_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta3 (degrees)');
    title('Theta3');
    grid on;
    
    %% V? Theta_dot
    axes(handles.plotTheta1Dot);
    theta1_dot_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta1 Dot (degrees/s)');
    title('Theta1 Dot');
    grid on;

    axes(handles.plotTheta2Dot);
    theta2_dot_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta2 Dot (degrees/s)');
    title('Theta2 Dot ');
    grid on;

    axes(handles.plotTheta3Dot);
    theta3_dot_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta3 Dot (degrees/s)');
    title('Theta3 Dot');
    grid on;
    %% Simulate Robot base on Velocity 
    for i = 1:length(time)
        [x, y, z] = forward_kinematics(theta1(i), theta2(i), theta3(i));
        trajectory_x(i) = x;
        trajectory_y(i) = y;
        trajectory_z(i) = z; 
        plot_robot_simulate(handles, theta1(i), theta2(i), theta3(i)); 
        view(60, 60);   
        theta1_deg = rad2deg(theta1(i)); 
        theta2_deg = rad2deg(theta2(i)); 
        theta3_deg = rad2deg(theta3(i)); 
        set(handles.theta1, 'String', sprintf('%.6f', theta1_deg));
        set(handles.theta2, 'String', sprintf('%.6f', theta2_deg));
        set(handles.theta3, 'String', sprintf('%.6f', theta3_deg));
       
        set(handles.edit_x, 'String', sprintf('%.6f', x));
        set(handles.edit_y, 'String', sprintf('%.6f', y));
        set(handles.edit_z, 'String', sprintf('%.6f', z));
        
        plot3(trajectory_x(1:i), trajectory_y(1:i), trajectory_z(1:i), 'r-', 'LineWidth', 2);
        
               



        
        
        
        %% Update data cho Theta
        set(theta1_plot, 'XData', time(1:i), 'YData', rad2deg(theta1(1:i)));  % ??
        set(theta2_plot, 'XData', time(1:i), 'YData', rad2deg(theta2(1:i)));  % ??
        set(theta3_plot, 'XData', time(1:i), 'YData', rad2deg(theta3(1:i)));  % ??
        set(theta1_plot, 'YData', rad2deg(theta1(1:i)));  % ??
        set(theta2_plot, 'YData', rad2deg(theta2(1:i)));  % ??
        set(theta3_plot, 'YData', rad2deg(theta3(1:i)));  % ??

        if i>1 
        %% Tính Theta_dot (V?n t?c góc)
        theta1_dot_deg = rad2deg((theta1(i) - theta1(i-1)) / (time(i) - time(i-1)));  % ??/giây
        theta2_dot_deg = rad2deg((theta2(i) - theta2(i-1)) / (time(i) - time(i-1)));  % ??/giây
        theta3_dot_deg = rad2deg((theta3(i) - theta3(i-1)) / (time(i) - time(i-1)));  % ??/giây
        %% L?u v?n t?c góc vào m?ng
        theta1_dot_array(i) = theta1_dot_deg;
        theta2_dot_array(i) = theta2_dot_deg;
        theta3_dot_array(i) = theta3_dot_deg;
        %% Update data cho Theta_dot
        set(theta1_dot_plot, 'XData', time(1:i), 'YData', theta1_dot_array(1:i));  % ??/giây
        set(theta2_dot_plot, 'XData', time(1:i), 'YData', theta2_dot_array(1:i));  % ??/giây
        set(theta3_dot_plot, 'XData', time(1:i), 'YData', theta3_dot_array(1:i));  % ??/giây
        set(theta1_dot_plot, 'YData', theta1_dot_array(1:i));  % ??/giây
        set(theta2_dot_plot, 'YData', theta2_dot_array(1:i));  % ??/giây
        set(theta3_dot_plot, 'YData', theta3_dot_array(1:i));  % ??/giây
        end
        drawnow;
        %% Check
        if i <= num_steps_acc
        fprintf('Step %d (Acceleration): Time = %.3f s, Velocity = %.3f , Acceleration = %.3f \n', ...
                i, time(i), velocity(i), a_max);
        elseif i <= num_steps_acc + num_steps_const
            fprintf('Step %d (Constant Velocity): Time = %.3f s, Velocity = %.3f \n', ...
                    i, time(i), velocity(i));
        else
            fprintf('Step %d (Deceleration): Time = %.3f s, Velocity = %.3f , Deceleration = %.3f \n', ...
                    i, time(i), velocity(i), -a_max);
        end
    end
   
axes(handles.axes1);  % Set current axes to axes1  % Clear the axes
hold on;  % Retain current plot to add new data
xlabel('X'); ylabel('Y'); zlabel('Z');  % Label the axes
grid on;  % Turn on the grid
plot3(trajectory_x, trajectory_y, trajectory_z, 'r-', 'LineWidth', 2);
hold off;

function x_a_Callback(hObject, eventdata, handles)
% hObject    handle to x_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_a as text
%        str2double(get(hObject,'String')) returns contents of x_a as a double


% --- Executes during object creation, after setting all properties.
function x_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_a_Callback(hObject, eventdata, handles)
% hObject    handle to y_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_a as text
%        str2double(get(hObject,'String')) returns contents of y_a as a double


% --- Executes during object creation, after setting all properties.
function y_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_a_Callback(hObject, eventdata, handles)
% hObject    handle to z_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_a as text
%        str2double(get(hObject,'String')) returns contents of z_a as a double


% --- Executes during object creation, after setting all properties.
function z_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_b_Callback(hObject, eventdata, handles)
% hObject    handle to x_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_b as text
%        str2double(get(hObject,'String')) returns contents of x_b as a double


% --- Executes during object creation, after setting all properties.
function x_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_b_Callback(hObject, eventdata, handles)
% hObject    handle to y_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_b as text
%        str2double(get(hObject,'String')) returns contents of y_b as a double


% --- Executes during object creation, after setting all properties.
function y_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_b_Callback(hObject, eventdata, handles)
% hObject    handle to z_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_b as text
%        str2double(get(hObject,'String')) returns contents of z_b as a double


% --- Executes during object creation, after setting all properties.
function z_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_m_Callback(hObject, eventdata, handles)
% hObject    handle to x_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_m as text
%        str2double(get(hObject,'String')) returns contents of x_m as a double


% --- Executes during object creation, after setting all properties.
function x_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_m_Callback(hObject, eventdata, handles)
% hObject    handle to y_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_m as text
%        str2double(get(hObject,'String')) returns contents of y_m as a double


% --- Executes during object creation, after setting all properties.
function y_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_m_Callback(hObject, eventdata, handles)
% hObject    handle to z_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_m as text
%        str2double(get(hObject,'String')) returns contents of z_m as a double


% --- Executes during object creation, after setting all properties.
function z_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_n_Callback(hObject, eventdata, handles)
% hObject    handle to x_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_n as text
%        str2double(get(hObject,'String')) returns contents of x_n as a double


% --- Executes during object creation, after setting all properties.
function x_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_n_Callback(hObject, eventdata, handles)
% hObject    handle to y_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_n as text
%        str2double(get(hObject,'String')) returns contents of y_n as a double


% --- Executes during object creation, after setting all properties.
function y_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_n_Callback(hObject, eventdata, handles)
% hObject    handle to z_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_n as text
%        str2double(get(hObject,'String')) returns contents of z_n as a double


% --- Executes during object creation, after setting all properties.
function z_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonCalc.
function buttonCalc_Callback(hObject, eventdata, handles)
% hObject    handle to buttonCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %point A
    x1 = str2double(get(handles.x_a,'String'));
    y1 = str2double(get(handles.y_a,'String'));
    z1 = str2double(get(handles.z_a,'String'));
    %point B
    x2 = str2double(get(handles.x_b,'String'));
    y2 = str2double(get(handles.y_b,'String'));
    z2 = str2double(get(handles.z_b,'String'));
    % Calculate distance
    d = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2); 
    % Parameter Vehicity max and ACC max
    v_max = 0.175; 
    a_max = 0.05;  
    % Calculate Time ACC and DEC 
    t_acc = sqrt(d / (2 * a_max));  %use formulas "d = 0.5 * a * t^2"
    t_dec = t_acc; 
    % Calculate distance AM:=> ACC ;  MN:=> CONST ; NB:=>DEC
    d_acc = 0.5 * a_max * t_acc^2; % d_AM
    d_dec = 0.5 * a_max * t_dec^2; % d_NB
    d_constant = d - d_acc - d_dec; % d_MN   
    %Calculate M; N 
     d_AM = 0.5 * a_max * t_acc^2;
     d_NB = d_AM;
    %point M
     x3 = x1 + (d_AM / d) * (x2 - x1);
     y3 = y1 + (d_AM / d) * (y2 - y1); 
     z3 = z1 + (d_AM / d) * (z2 - z1); 
    %point N
     x4 = x1 + ((d_AM + d_constant) / d) * (x2 - x1);
     y4 = y1 + ((d_AM + d_constant) / d) * (y2 - y1);
     z4 = z1 + ((d_AM + d_constant) / d) * (z2 - z1);
    
   set(handles.x_m, 'String', sprintf('%.6f', x3));
   set(handles.y_m, 'String', sprintf('%.6f', y3));
   set(handles.z_m, 'String', sprintf('%.6f', z3));
   set(handles.x_n, 'String', sprintf('%.6f', x4));
   set(handles.y_n, 'String', sprintf('%.6f', y4));
   set(handles.z_n, 'String', sprintf('%.6f', z4));
% --- Executes on button press in buttonMove.
function buttonMove_Callback(hObject, eventdata, handles)
   %point A
    x1 = str2double(get(handles.x_a,'String'));
    y1 = str2double(get(handles.y_a,'String'));
    z1 = str2double(get(handles.z_a,'String'));
    %point B
    x2 = str2double(get(handles.x_b,'String'));
    y2 = str2double(get(handles.y_b,'String'));
    z2 = str2double(get(handles.z_b,'String'));
    % Calculate distance
    d = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2); 
    % Parameter Vehicity max and ACC max
    
    v_max = str2double(get(handles.edit_v,'String')); 
    a_max = str2double(get(handles.edit_a,'String'));  
    % Calculate Time ACC and DEC 
    t_acc = sqrt(d / (2 * a_max));  %use formulas "d = 0.5 * a * t^2"
    t_dec = t_acc; 
    % Calculate distance AM:=> ACC ;  MN:=> CONST ; NB:=>DEC
    d_acc = 0.5 * a_max * t_acc^2; % d_AM
    d_dec = 0.5 * a_max * t_dec^2; % d_NB
    d_constant = d - d_acc - d_dec; % d_MN   
    %Calculate M; N 
     d_AM = 0.5 * a_max * t_acc^2;
     d_NB = d_AM;
    %point M
     x3 = x1 + (d_AM / d) * (x2 - x1);
     y3 = y1 + (d_AM / d) * (y2 - y1); 
     z3 = z1 + (d_AM / d) * (z2 - z1); 
    %point N
     x4 = x1 + ((d_AM + d_constant) / d) * (x2 - x1);
     y4 = y1 + ((d_AM + d_constant) / d) * (y2 - y1);
     z4 = z1 + ((d_AM + d_constant) / d) * (z2 - z1);
    % Display point M, N position
%    set(handles.x_m, 'String', sprintf('%.6f', x3));
%    set(handles.y_m, 'String', sprintf('%.6f', y3));
%    set(handles.z_m, 'String', sprintf('%.6f', z3));
%    set(handles.x_n, 'String', sprintf('%.6f', x4));
%    set(handles.y_n, 'String', sprintf('%.6f', y4));
%    set(handles.z_n, 'String', sprintf('%.6f', z4));
%   use Inverse kinematics calculate theta
    [thetaA(1), thetaA(2), thetaA(3)] = inverse_kinematics(x1,y1,z1);
    [thetaB(1), thetaB(2), thetaB(3)] = inverse_kinematics(x2,y2,z2);
    [thetaM(1), thetaM(2), thetaM(3)] = inverse_kinematics(x3,y3,z3);
    [thetaN(1), thetaN(2), thetaN(3)] = inverse_kinematics(x4,y4,z4);
    
    
    
    if d_constant > 0
        t_constant = d_constant / v_max; 
    else
        t_constant = 0;
    end
    % Total time from A to B
    t_total = t_acc + t_constant + t_dec;
    % step for draw
    num_steps = 100; 
    time_acc = linspace(0, t_acc, num_steps * (t_acc / t_total));  
    time_const = linspace(t_acc, t_acc + t_constant, num_steps * (t_constant / t_total)); 
    time_dec = linspace(t_acc + t_constant, t_total, num_steps * (t_dec / t_total));  
   
    v_acc = a_max * time_acc; 
    v_const = v_max * ones(1, length(time_const));  
    v_dec = v_max - a_max * (time_dec - (t_acc + t_constant));  
    %% Draw Velocity of Robot
    time = [time_acc, time_const, time_dec];  
    velocity = [v_acc, v_const, v_dec];  
    axes(handles.velocity);  
    plot(time, velocity, 'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('Velocity');
    title('Velocity');
    grid on;
    %% Draw Position 
    pos_acc = 0.5 * a_max * time_acc.^2; 
    pos_const = d_acc + v_max * (time_const - t_acc);  
    pos_dec = d_acc + d_constant + v_max * (time_dec - (t_acc + t_constant)) - 0.5 * a_max * (time_dec - (t_acc + t_constant)).^2; 
    position = [pos_acc, pos_const, pos_dec];
    axes(handles.position);  
    plot(time, position, 'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('Position');
    title('Position');
    grid on;
    %% Draw Acceleration
    acceleration_acc = a_max * ones(1, length(time_acc));
    acceleration_const = zeros(1, length(time_const));
    acceleration_dec =  -a_max * ones(1, length(time_dec));
    acceleration = [acceleration_acc, acceleration_const, acceleration_dec];
    
    axes(handles.acceleration);  
    plot(time, acceleration, 'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('acceleration');
    title('Acceleration');
    grid on;
    %% Time steps
    num_steps_acc = round(num_steps * (t_acc / t_total));        
    num_steps_const = round(num_steps * (t_constant / t_total));  
    num_steps_dec = num_steps - num_steps_acc - num_steps_const; 
    %% ACC 
    theta1_acc = linspace(thetaA(1), thetaM(1), num_steps_acc);
    theta2_acc = linspace(thetaA(2), thetaM(2), num_steps_acc);
    theta3_acc = linspace(thetaA(3), thetaM(3), num_steps_acc);
    %% CONST
    theta1_const = linspace(thetaM(1), thetaN(1), num_steps_const);  
    theta2_const = linspace(thetaM(2), thetaN(2), num_steps_const);
    theta3_const = linspace(thetaM(3), thetaN(3), num_steps_const);
    %% DEC
    theta1_dec = linspace(thetaN(1), thetaB(1), num_steps_dec); 
    theta2_dec = linspace(thetaN(2), thetaB(2), num_steps_dec);
    theta3_dec = linspace(thetaN(3), thetaB(3), num_steps_dec);
    %% Final
    theta1 = [theta1_acc, theta1_const, theta1_dec];
    theta2 = [theta2_acc, theta2_const, theta2_dec];
    theta3 = [theta3_acc, theta3_const, theta3_dec];

    %% Simulate Robot base on Velocity 
    trajectory_x = zeros(1, length(time));
    trajectory_y = zeros(1, length(time));
    trajectory_z = zeros(1, length(time));
    

    %% Kh?i t?o m?ng
    theta1_dot_array = zeros(1, length(time));
    theta2_dot_array = zeros(1, length(time));
    theta3_dot_array = zeros(1, length(time));
        
    %% V? Theta
    axes(handles.plotTheta1);
    theta1_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta1 (degrees)');
    title('Theta1');
    grid on;

    axes(handles.plotTheta2);
    theta2_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta2 (degrees)');
    title('Theta2');
    grid on;
  
    axes(handles.plotTheta3);
    theta3_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta3 (degrees)');
    title('Theta3');
    grid on;
    %% V? Theta_dot
    axes(handles.plotTheta1Dot);
    theta1_dot_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta1 Dot (degrees/s)');
    title('Theta1 Dot');
    grid on;

    axes(handles.plotTheta2Dot);
    theta2_dot_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta2 Dot (degrees/s)');
    title('Theta2 Dot');
    grid on;

    axes(handles.plotTheta3Dot);
    theta3_dot_plot = plot(time, zeros(1, length(time)), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Theta3 Dot (degrees/s)');
    title('Theta3 Dot');
    grid on;

%% Vòng lap cho simulate: 
    for i = 1:length(time)
        [x, y, z] = forward_kinematics(theta1(i), theta2(i), theta3(i));
        plot_robot_simulate(handles, theta1(i), theta2(i), theta3(i)); 
        view(60, 60); 
        theta1_deg = rad2deg(theta1(i)); 
        theta2_deg = rad2deg(theta2(i)); 
        theta3_deg = rad2deg(theta3(i)); 
        %% Update data cho Theta
        set(theta1_plot, 'XData', time(1:i), 'YData', rad2deg(theta1(1:i)));  
        set(theta2_plot, 'XData', time(1:i), 'YData', rad2deg(theta2(1:i)));  
        set(theta3_plot, 'XData', time(1:i), 'YData', rad2deg(theta3(1:i)));  
        set(theta1_plot, 'YData', rad2deg(theta1(1:i)));  
        set(theta2_plot, 'YData', rad2deg(theta2(1:i)));  
        set(theta3_plot, 'YData', rad2deg(theta3(1:i)));  
        if i > 1
        %% Tính Theta_dot (van toc góc)
        theta1_dot_deg = rad2deg((theta1(i) - theta1(i-1)) / (time(i) - time(i-1)));  % ??/giây
        theta2_dot_deg = rad2deg((theta2(i) - theta2(i-1)) / (time(i) - time(i-1)));  % ??/giây
        theta3_dot_deg = rad2deg((theta3(i) - theta3(i-1)) / (time(i) - time(i-1)));  % ??/giây
        %% Luu van toc góc vào mang
        theta1_dot_array(i) = theta1_dot_deg;
        theta2_dot_array(i) = theta2_dot_deg;
        theta3_dot_array(i) = theta3_dot_deg;
        %% Update data cho Theta_dot
        set(theta1_dot_plot, 'XData', time(1:i), 'YData', theta1_dot_array(1:i));  % ??/giây
        set(theta2_dot_plot, 'XData', time(1:i), 'YData', theta2_dot_array(1:i));  % ??/giây
        set(theta3_dot_plot, 'XData', time(1:i), 'YData', theta3_dot_array(1:i));  % ??/giây
        set(theta1_dot_plot, 'YData', theta1_dot_array(1:i));  % ??/giây
        set(theta2_dot_plot, 'YData', theta2_dot_array(1:i));  % ??/giây
        set(theta3_dot_plot, 'YData', theta3_dot_array(1:i));  % ??/giây
        end
        set(handles.theta1, 'String', sprintf('%.6f', theta1_deg));
        set(handles.theta2, 'String', sprintf('%.6f', theta2_deg));
        set(handles.theta3, 'String', sprintf('%.6f', theta3_deg));    
        set(handles.edit_x, 'String', sprintf('%.6f', x));
        set(handles.edit_y, 'String', sprintf('%.6f', y));
        set(handles.edit_z, 'String', sprintf('%.6f', z));
        
        trajectory_x(i) = x;
        trajectory_y(i) = y;
        trajectory_z(i) = z;
        plot3(trajectory_x(1:i), trajectory_y(1:i), trajectory_z(1:i), 'r-', 'LineWidth', 2);
        drawnow;
        %% Check 
        if i <= num_steps_acc
        fprintf('Step %d (Acceleration): Time = %.3f s, Velocity = %.3f , Acceleration = %.3f \n', ...
                i, time(i), velocity(i), a_max);
        elseif i <= num_steps_acc + num_steps_const
            fprintf('Step %d (Constant Velocity): Time = %.3f s, Velocity = %.3f \n', ...
                    i, time(i), velocity(i));
        else
            fprintf('Step %d (Deceleration): Time = %.3f s, Velocity = %.3f , Deceleration = %.3f \n', ...
                    i, time(i), velocity(i), -a_max);
        end
    end
    %% Draw Line
    axes(handles.axes1);  % Set current axes to axes1  % Clear the axes
    hold on;  % Retain current plot to add new data
    xlabel('X'); ylabel('Y'); zlabel('Z');  % Label the axes
    grid on;  % Turn on the grid
    plot3(trajectory_x, trajectory_y, trajectory_z, 'r-', 'LineWidth', 2);    
    hold off;
% --- Executes on button press in plotWorkspaceCheckbox.
function plotWorkspaceCheckbox_Callback(hObject, eventdata, handles)
     plotWorkspace = get(hObject, 'Value');
     handles.plotWorkspace = plotWorkspace;
     guidata(hObject, handles);
% --- Executes on button press in plotCoordinatesCheckbox.
function plotCoordinatesCheckbox_Callback(hObject, eventdata, handles)
    plotCoordinate = get(hObject, 'Value');
    handles.plotCoordinate = plotCoordinate;
    guidata(hObject, handles);
function edit_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_v as text
%        str2double(get(hObject,'String')) returns contents of edit_v as a double
% --- Executes during object creation, after setting all properties.
function edit_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_a as text
%        str2double(get(hObject,'String')) returns contents of edit_a as a double
% --- Executes during object creation, after setting all properties.
function edit_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
