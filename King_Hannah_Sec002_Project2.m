%% EMEC303
% Project 2
% Team: Trajectory
% Projectile Motion and Trajectory
% Hannah King

clear; clc;
%% Inputs
% Got skeleton code for this part from Mathworks

%Initialize Variables
xloc = 0; % x location of target m
yloc = 0; % y location of target m
zloc = 0; % z location of target m
v = 0; % initial velocity of projectile m/s

% Ask user for initial values.
defaultValue = {num2str(xloc), num2str(yloc),...
    num2str(zloc), num2str(v)};
titleBar = 'Enter values';
userPrompt = {'Location of Target (Front To Back) in meters : ',... 
    'Location of Target (Height) in meters : ',... 
    'Location of Target (Side to Side) in meters : ',...
    'Initial Velocity in meters per second : '};
UserInput = inputdlg(userPrompt, titleBar, 1,...
    defaultValue, 'on');

% Cancel Button
if isempty(UserInput),return,end 

% Convert to floating point from string.
usersValue1 = str2double(UserInput{1});
usersValue2 = str2double(UserInput{2});
usersValue3 = str2double(UserInput{3});
usersValue4 = str2double(UserInput{4});

% Check for a valid number.
if isnan(usersValue1)
    usersValue1 = str2double(defaultValue{1});
    message = sprintf('This is not a valid number, please try again.');
    uiwait(warndlg(message));
    UserInput = inputdlg(userPrompt, titleBar, 1, defaultValue, 'on');
end

if isnan(usersValue2)
    usersValue2 = str2double(defaultValue{2});
    message = sprintf('This is not a valid number, please try again.');
    uiwait(warndlg(message));
    UserInput = inputdlg(userPrompt, titleBar, 1, defaultValue, 'on');
end

if isnan(usersValue3)
    usersValue3 = str2double(defaultValue{3});
    message = sprintf('This is not a valid number, please try again.');
    uiwait(warndlg(message));
    UserInput = inputdlg(userPrompt, titleBar, 1, defaultValue, 'on');
end

% Assign user's values:
xloc = usersValue1; % m
yloc = usersValue2; % m
zloc= usersValue3; % m
v = usersValue4; % m/s

%% Computations
%Variables 
g = 9.81; % Gravity m/s^2
Cd = 0.507; % Coef of Drag for Tennis ball
Cl = 0.3; % Coeff of Lift for Tennis ball (Backspin)
p = 1.225; % Density of air kg/m
A = 0.0034; % Cross secional area of Tennis Ball m^2
m = 0.0577; % Mass of Tennis Ball kg
angle_min = 15; % Minimum angle robot can shoot
angle_max = 75; % Maximum angle robot can shoot
s = 0.45; % Angle step size (1/4 motor revolutiondue to 4:1 gear ratio)
h = 0.01; % Time step
k = 0.01; % Velocity step
xi = 0; % Initial x location m
yi = 0.33118; % Initial y location m
r = 0.1524/2; % Target radius m

d = sqrt(xloc^2 + zloc^2); % Distance from Robot to Target
hangle = atan(xloc/zloc); % Angle of rotation in the horizontal plane

L = 0.5 * p * Cl * A; % Lift Force N
D = 0.5 * p * Cd * A; % Drag Force N
G = g * m; %Gravitational Force N

%Iterate
vangle = angle_min; %Initial angle of rotation in the vertical plane
Check = 0;
hit = 0;

while Check == 0 % Target not hit
    xf = 0; % Initialize x_final location
    yf = 0; % Initialize x_final location
    t = 0; % Initialize time s
    xa = []; % Initialize x-adjusted array
    y = []; % Initialize y-loc array
    
    Vx = v*cosd(vangle); % Velocity in x direction
    Vy = v*sind(vangle); % Velocity in y direction

    Lx = (L * v^2) * sind(vangle); % Lift in the x direction
    Ly = (L * v^2) * cosd(vangle); % Lift in the y direction
    Dx = (D * v^2) * cosd(vangle); % Drag in the x direction
    Dy = (D * v^2) * sind(vangle); % Drag in the y direction

    ax = (-Dx - Lx)/m; % Accel in the x direction m/s^2
    ay = (-Dy + Ly - G)/m; % Accel in the y direction m/s^2
    
    % Position Calculation
    while xf < (d - r) % The ball hasn't reached the target
        if xf >= 0 % The ball's location is positive
            xf = xi + (Vx * t) + (0.5 * ax * t^2); % Calculating x-loc     
            xa = [xa,xf]; % Appending the x-location array 
            yf = yi + (Vy * t) + (0.5 * ay * t^2); % Calculating y-loc
            y = [y,yf]; % Appending the y-;ocation array
            t = t + h; % Time update
        else
            break
        end
    end
    
    % Angle adjustment
    if yf > (yloc+r) || yf < (yloc-r) % Target not hit in the y direction
        if vangle < (angle_max - s) % Angle is under max angle
            vangle = vangle + s; % Increase angle
        else % Target cannot be hit under these conditions
            Check = 1;
            disp('Target cannot be hit with the current parameters.');
            disp('Shooting at a 42 degree angle, a new initial velocity will be calculated.')
            disp(' ')
        end   
    else % Target hit
        Check = 1;
        disp("Target Hit!")
        hit = 1;
        tf = t;
    end
end

if hit == 1 % Print for Initial hit
    print1 = ('The horizontal angle is %4.2f degrees.');
    fprintf(print1, hangle); 
    disp(' ');
    print2 = ('The vertical angle is %4.2f degrees.');
    fprintf(print2, vangle);
    disp(' ');
    print3 = ('The time is takes for the Tennis Ball to hit the Target is %4.2f seconds');
    fprintf(print3, tf);
    
else % Velocity Change Required
    vangle = 42; % Default angle
    Vo = v; % Set initial velocity
    
    if yf > (yloc+r) % Overshooting target
        while yf > (yloc+r) % Loop
            xf = 0; % Initialize x_final location
            yf = 0; % Initialize x_final location
            t = 0; % initiaize time

            xa = []; % Initialize x-adjusted array
            y = []; % Initialize y-loc array
            
            if Vo > 0
                Vo = Vo - k; % Decrease Velocity
            else
                disp("Cannot reach the target.")
            end
            
            Vx = Vo*cosd(vangle); % Velocity in x direction
            Vy = Vo*sind(vangle); % Velocity in y direction

            Lx = (L * Vo^2) * sind(vangle); % Lift in the x direction
            Ly = (L * Vo^2) * cosd(vangle); % Lift in the y direction
            Dx = (D * Vo^2) * cosd(vangle); % Drag in the x direction
            Dy = (D * Vo^2) * sind(vangle); % Drag in the y direction

            ax = (-Dx - Lx)/m; % Accel in the x direction m/s^2
            ay = (-Dy + Ly - G)/m; % Accel in the y direction m/s^2
       
            % Check Positions
            while xf < (d - r) % The ball hasn't reached the target
                if xf >= 0 % The ball's location is positive
                    xf = xi + (Vx * t) + (0.5 * ax * t^2); % Calculating x-loc     
                    xa = [xa,xf]; % Appending the x-location array 
                    yf = yi + (Vy * t) + (0.5 * ay * t^2); % Calculating y-loc
                    y = [y,yf]; % Appending the y-;ocation array
                    t = t + h; % Time update
                else
                    break
                end
                tf = t;
            end
        end
        
        % Print for adjusted hit
        disp("Target Hit with adjustments.")
        print1 = ('The horizontal angle is %4.2f degrees.');
        fprintf(print1, hangle); 
        disp(' ');
        print2 = ("The initial velocity is %4.2f meters per second.");
        fprintf(print2, Vo);
        disp(' ');
        print3 = ('The time is takes for the Tennis Ball to hit the Target is %4.2f seconds');
        fprintf(print3, tf);
        
    elseif yf < (yloc-r) % Undershooting
        while yf < (yloc-r) % Loop
            xf = 0; % Initialize x_final location
            yf = 0; % Initialize x_final location
            t = 0; % Initialize time

            xa = []; % Initialize x-adjusted array
            y = []; % Initialize y-loc array
    
            Vo = Vo + k; % Increase Velocity
            
            Vx = Vo*cosd(vangle); % Velocity in x direction
            Vy = Vo*sind(vangle); % Velocity in y direction

            Lx = (L * Vo^2) * sind(vangle); % Lift in the x direction
            Ly = (L * Vo^2) * cosd(vangle); % Lift in the y direction
            Dx = (D * Vo^2) * cosd(vangle); % Drag in the x direction
            Dy = (D * Vo^2) * sind(vangle); % Drag in the y direction

            ax = (-Dx - Lx)/m; % Accel in the x direction m/s^2
            ay = (-Dy + Ly - G)/m; % Accel in the y direction m/s^2
       
            %Check
            while xf < (d - r) % The ball hasn't reached the target
                if xf >= 0 % The ball's location is positive
                    xf = xi + (Vx * t) + (0.5 * ax * t^2); % Calculating x-loc     
                    xa = [xa,xf]; % Appending the x-location array 
                    yf = yi + (Vy * t) + (0.5 * ay * t^2); % Calculating y-loc
                    y = [y,yf]; % Appending the y-;ocation array
                    t = t + h; % Time update
                else
                    break
                end
            end
            tf = t;
        end
        
        % Print for adjusted hit
        disp("Target hit with adjustments.")
        print1 = ('The horizontal angle is %4.2f degrees.');
        fprintf(print1, hangle); 
        disp(' ');
        print2 = ("The initial velocity is %4.2f meters per second.");
        fprintf(print2, Vo);
        disp(' ');
        print3 = ('The time is takes for the Tennis Ball to hit the Target is %4.2f seconds');
        fprintf(print3, tf);
        
    end
end

%% Plot Trajectory
% 2D plot
figure(1); clf(1);
plot(xa, y)
hold on
title("Tennis Ball Trajectory")
xlabel("Distance (m)")
ylabel("Height (m)")
grid on
plot(d, yloc, 'ro', 'MarkerSize', 30)

% 3D plot
x = []; % Initialize x-loc array
z = []; % Initialize z-loc array
i = 1;

for i=1:length(xa) % Break xa up into x and z components
    xtrue = xa(i) * sin(hangle);
    x = [x, xtrue];
    ztrue = xa(i) * cos(hangle);
    z = [z, ztrue];
    i = i + 1;
end

figure(2); clf(2);
plot3(x, z, y)
hold on
title("Tennis Ball Trajectory")
xlabel("Distance: Front-to-Back (m)")
ylabel('Distance: Side-to-Side (m)')
zlabel("Height (m)")
grid on
plot3(xloc, zloc, yloc, 'ro', 'MarkerSize', 30)


