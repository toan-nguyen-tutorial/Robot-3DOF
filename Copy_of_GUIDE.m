function varargout = GUIDE(varargin)
% GUIDE MATLAB code for GUIDE.fig
%      GUIDE, by itself, creates a new GUIDE or raises the existing
%      singleton*.
%
%      H = GUIDE returns the handle to a new GUIDE or the handle to
%      the existing singleton*.
%
%      GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE.M with the given input arguments.
%
%      GUIDE('Property','Value',...) creates a new GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE

% Last Modified by GUIDE v2.5 23-Nov-2024 14:37:05

% Begin initialization code - DO NOT EDIT
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
% End initialization code - DO NOT EDIT


% --- Executes just before GUIDE is made visible.
function GUIDE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIDE (see VARARGIN)
% Setup the initial interface
axes(handles.axes1);
cla;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([-2 2 -2 2 -0.5 2]);
view(3);

% Set initial values
set(handles.theta1, 'String', '0');
set(handles.theta2, 'String', '0');
set(handles.theta3, 'String', '0');
set(handles.edit_x, 'String', '0.00');
set(handles.edit_y, 'String', '0.00');
set(handles.edit_z, 'String', '0.00');

% Plot robot at the default position
plot_robot(handles, 0, 0, 0);

% Choose default command line output for GUIDE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIDE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function theta1_Callback(hObject, eventdata, handles)
% hObject    handle to theta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta1 as text
%        str2double(get(hObject,'String')) returns contents of theta1 as a double


% --- Executes during object creation, after setting all properties.
function theta1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta2_Callback(hObject, eventdata, handles)
% hObject    handle to theta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta2 as text
%        str2double(get(hObject,'String')) returns contents of theta2 as a double


% --- Executes during object creation, after setting all properties.
function theta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta3_Callback(hObject, eventdata, handles)
% hObject    handle to theta3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta3 as text
%        str2double(get(hObject,'String')) returns contents of theta3 as a double


% --- Executes during object creation, after setting all properties.
function theta3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_x as text
%        str2double(get(hObject,'String')) returns contents of edit_x as a double


% --- Executes during object creation, after setting all properties.
function edit_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_y as text
%        str2double(get(hObject,'String')) returns contents of edit_y as a double


% --- Executes during object creation, after setting all properties.
function edit_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_z as text
%        str2double(get(hObject,'String')) returns contents of edit_z as a double


% --- Executes during object creation, after setting all properties.
function edit_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in button_fwd.
function button_fwd_Callback(hObject, eventdata, handles)
% hObject    handle to button_fwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get input angles
theta1 = str2double(get(handles.theta1, 'String'));
theta2 = str2double(get(handles.theta2, 'String'));
theta3 = str2double(get(handles.theta3, 'String'));

if isnan(theta1) || isnan(theta2) || isnan(theta3)
    msgbox('Please enter valid angles for Theta1, Theta2, and Theta3.', 'Error', 'error');
    return;
end

% Convert to radians
theta1_rad = deg2rad(theta1);
theta2_rad = deg2rad(theta2);
theta3_rad = deg2rad(theta3);
% Forward kinematics to calculate X, Y, Z
[x, y, z] = forward_kinematics(theta1_rad, theta2_rad, theta3_rad);

% Display results
set(handles.edit_x, 'String', sprintf('%.4f', x));
set(handles.edit_y, 'String', sprintf('%.4f', y));
set(handles.edit_z, 'String', sprintf('%.4f', z));

% Update robot plot
plot_robot(handles, theta1_rad, theta2_rad, theta3_rad);


%% --- Executes on button press in button_inv.
function button_inv_Callback(hObject, eventdata, handles)
% hObject    handle to button_inv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get input coordinates
x = str2double(get(handles.edit_x, 'String'));
y = str2double(get(handles.edit_y, 'String'));
z = str2double(get(handles.edit_z, 'String'));

if isnan(x) || isnan(y) || isnan(z)
    msgbox('Please enter valid coordinates for X, Y, and Z.', 'Error', 'error');
    return;
end

try
   % Inverse kinematics to calculate Theta1, Theta2, Theta3
   [theta1, theta2, theta3] = inverse_kinematics(x, y, z);

   % Convert to degrees if necessary
   theta1_deg = rad2deg(theta1); % Convert radians to degrees
   theta2_deg = rad2deg(theta2);
   theta3_deg = rad2deg(theta3);

   % Display results
   set(handles.theta1, 'String', sprintf('%.4f', theta1_deg));
   set(handles.theta2, 'String', sprintf('%.4f', theta2_deg));
   set(handles.theta3, 'String', sprintf('%.4f', theta3_deg));

   % Update robot plot
   plot_robot(handles, theta1, theta2, theta3);
catch
   msgbox('Position is not reachable.', 'Error', 'error');
end

%% --- Forward kinematics function
function [x, y, z] = forward_kinematics(theta1, theta2, theta3)
a1 = 1; a2 = 1; d1 = 0.4;  % Lengths of links and offsets

    % Compute transformation matrices for each joint using DH convention
    T1 = DHMatrix(theta1, d1, 0, pi/2);  % Transformation from base to joint 1
    T2 = T1 * DHMatrix(theta2, 0, a1, 0);  % Transformation from joint 1 to joint 2
    T3 = T2 * DHMatrix(theta3, 0, a2, 0);  % Transformation from joint 2 to joint 3

    % Final position of the end-effector (P) in homogeneous coordinates
    P = T3 * [0; 0; 0; 1];  % Multiply with [0; 0; 0; 1] to get position

    % Extract the x, y, z components from the resulting vector
    x = P(1); 
    y = P(2); 
    z = P(3);
%% --- Inverse kinematics function
function [theta1, theta2, theta3] = inverse_kinematics(x, y, z)
    % Lengths of links and offsets
    a1 = 1;    % Length of the first link
    a2 = 1;    % Length of the second link
    a3 = 1;    % Length of the third link
    d1 = 0.4;  % Offset along the Z axis

    % Compute distance r in the xy plane
    r = sqrt(x^2 + y^2); 

    % Compute the vertical distance s in z direction
    s = z - d1;

    % Solve for D using the inverse kinematics equations
    D = (r^2 + s^2 - a1^2 - a2^2) / (2 * a1 * a2);

    % Check if the position is reachable
    if abs(D) > 1
        error('Position not reachable.');  % If D is out of bounds, it's unreachable
    end

    % Calculate the angles using inverse kinematics


    theta1 = atan2(y, x); 
    theta3 = atan2(sqrt(1 - D^2), D);
    c3 = cos(theta3);
    s3 = sin(theta3);

    % Calculate c2 and s2 for theta2
    c2 = (r * (a2 + a3 * c3) + s * a3 * s3) / (a2^2 + a3^2 + 2 * a2 * a3 * c3);
    s2 = (-r * a3 * s3 - s * (a2 + a3 * c3)) / (a2^2 + a3^2 + 2 * a2 * a3 * c3);

    % Calculate theta2 using atan2
    theta2 = atan2(s2, c2);
%% --- Denavit-Hartenberg matrix function
function T = DHMatrix(theta, d, a, alpha)
T = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
     0,           sin(alpha),             cos(alpha),            d;
     0,           0,                      0,                     1];

%% --- Robot plotting function
function plot_robot(handles, theta1, theta2, theta3)
a1 = 1; a2 = 1; d1 = 0.4;
T1 = DHMatrix(theta1, d1, 0, pi/2);
T2 = T1 * DHMatrix(theta2, 0, a1, 0);
T3 = T2 * DHMatrix(theta3, 0, a2, 0);
P0 = [0; 0; 0; 1];
P1 = T1 * P0;
P2 = T2 * P0;
P3 = T3 * P0;

axes(handles.axes1);

cla;
hold on;
grid on;
%% base 
% V? base robot (hình tr?)
    radius_base = 0.5;  % Bán kính base
    height_base = 0.2;  % Chi?u cao base
    color_base = [0.7 0.7 0.7];  % Màu s?c base
    VeHinhTru(handles, 0, 0, -height_base/2, radius_base, height_base, color_base);  % V? hình tr? làm base
% Plot links
plot3([P0(1) P1(1)], [P0(2) P1(2)], [P0(3) P1(3)], 'r', 'LineWidth', 2);
plot3([P1(1) P2(1)], [P1(2) P2(2)], [P1(3) P2(3)], 'g', 'LineWidth', 2);
plot3([P2(1) P3(1)], [P2(2) P3(2)], [P2(3) P3(3)], 'b', 'LineWidth', 2);

% Plot joints
plot3(P1(1), P1(2), P1(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot3(P2(1), P2(2), P2(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(P3(1), P3(2), P3(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

xlabel('X'); ylabel('Y'); zlabel('Z');
axis([-2 2 -2 2 -0.5 2]);
view(3);
rotate3d(handles.axes1, 'on');

hold off;



%% --- Function to draw a cuboid (rectangular box)
function VeHop(handles, x0, y0, z0, a, b, h, colr)
    % Draw a cuboid (rectangular box) with given dimensions and position
    % x0, y0, z0: base position coordinates
    % a, b, h: dimensions of the cuboid (length, width, height)
    % colr: color for the cuboid
    
    [X, Y] = meshgrid([-a/2, a/2], [-b/2, b/2]);
    Z = [z0*ones(size(X)); (z0 + h)*ones(size(X))];
    
    surf(handles, X, Y, Z, 'FaceColor', colr, 'EdgeColor', 'none');
    fill3(handles, X(1,:), Y(1,:), Z(1,:), colr);  % Bottom face
    fill3(handles, X(2,:), Y(2,:), Z(2,:), colr);  % Top face
%% --- Function to draw a cylinder
function VeHinhTru(handles, x0, y0, z0, r, h, colr)
    % Draw a cylinder with given radius, height, and base position
    % x0, y0, z0: position of the base center
    % r: radius of the cylinder
    % h: height of the cylinder
    % colr: color for the cylinder
    
    [X, Y, Z] = cylinder(r, 100);  % Create the cylinder
    X = X + x0;  % Translate the cylinder
    Y = Y + y0;
    Z = Z * h + z0;  % Adjust the height
    
    surf(handles, X, Y, Z, 'FaceColor', colr, 'EdgeColor', 'none');
    fill3(handles, X(1,:), Y(1,:), Z(1,:), colr);  % Bottom face of the cylinder
    fill3(handles, X(2,:), Y(2,:), Z(2,:), colr);  % Top face of the cylinder