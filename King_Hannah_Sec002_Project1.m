%% EMEC 303
%Project 1
%Hannah King

clear; clc;
%%
%Inputs
k1 = 0.2; %W/m^2K, Coeff Thermal Res room to outside
k3 = 0.5; %W/m^2K, Coeff Thermal Res roof to outside
A_in = 0.5; %m^2, Area of fireplace
C_air = 1005; %J/kgK
ro_air = 1.2; %kg/m^3
h = 1; %step size 1 sec
L = 5; %m
W = 6; %m
H1 = 3; %m
H2 = 3; %m
A1 = 2 * (H1*L) + 2 * (H1*W); %m^2, Surface area of the room walls
A2 = L * W; %m^2, Surface area between room and roof
A3 = (H2*W) + 2 * (4.24*L); %m^2, Surface area of the roof walls
V_room = (L * W * H1); %m^3, Volume of room
V_roof = ((H2 * W)/2) * L; %m^3, Voume of roof
t = [0:h:(86400*2)]; %s

%Variable Inputs
Q_in_low = 1250; %W/m^2, Heat flux from freplace
Q_in_high = 2450; %W/m^2, Heat flux from freplace
k2_closed = 2.5; %W/m^2K, Coeff Thermal Res room to roof
k2_open = 5; %W/m^2K, Coeff Thermal Res room to roof open door

%Functions
To = @(t) -10*sin(2*pi*t / 86400); %Toutside, deg C

%Initialize Arrays
T1 = zeros(1, length(t)); %deg C, Temp room
T2 = zeros(1, length(t)); %deg C, Temp roof
Q1 = zeros(1, length(t)); %W/m^2, Heat transfer from room to outside
Q2 = zeros(1, length(t)); %W/m^2, Heat transfer from room to roof
Q3 = zeros(1, length(t)); %W/m^2, Heat transfer from roof to outside

%Initial conditions
T1(1) = 5; %deg C, room
T2(1) = 7; %deg C, roof

%Iterate
for i = 1:length(t)-1
    %time
    s = t(i);
    
    %Heat fluxes
    Q1(i) = k1 * (T1(i) - To(t(i)));
    Q3(i) = k3 * (T2(i) - To(t(i)));
    
    if s < 40000 
        %Heat flux
        Q2(i) = k2_open * (T1(i) - T2(i) + 5); 
        
        %dT/dt equations
        deltaT1 = ((A_in*Q_in_high) - (A1*Q1(i)) - (A2*Q2(i))) /...
            (C_air * ro_air * V_room); %dec C, Temp of room
        deltaT2 = ((A2*Q2(i)) - (A3*Q3(i))) /...
            (C_air * ro_air * V_roof); %dec C, Temp of roof
    
    elseif 87500 < s && s < 127500
        %Heat flux
        Q2(i) = k2_open * (T1(i) - T2(i) + 5); 
        
        %dT/dt equations
        deltaT1 = ((A_in*Q_in_high) - (A1*Q1(i)) - (A2*Q2(i))) /...
            (C_air * ro_air * V_room); %dec C, Temp of room
        deltaT2 = ((A2*Q2(i)) - (A3*Q3(i))) /...
            (C_air * ro_air * V_roof); %dec C, Temp of roof
    else
        %Heat flux
        Q2(i) = k2_closed * (T1(i) - T2(i) + 5); 
        
        %dT/dt equations
            deltaT1 = ((A_in*Q_in_low) - (A1*Q1(i)) - (A2*Q2(i))) /...
                (C_air * ro_air * V_room); %dec C, Temp of room
            deltaT2 = ((A2*Q2(i)) - (A3*Q3(i))) /...
                (C_air * ro_air * V_roof); %dec C, Temp of roof
    end
    
    %Euler's
    T1(i+1) = T1(i) + h*deltaT1;
    T2(i+1) = T2(i) + h*deltaT2;
   
end

figure(1)
clf(1)
plot(t, T1)
hold on 
plot(t, T2)
plot(t, To(t))
xlabel('Time (sec)')
ylabel('Temperature (deg C)')
title('Temperature Over Time')
legend('T1', 'T2', 'To')