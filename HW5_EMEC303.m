%% EMEC303
% HW 5
% Hannah King

clear; clc;
%% Problem 1
%Inputs
m = 5; %kg
k = 0.5; %N/m
wo = sqrt(k/m); 
A = 1;
B = 1;

%Initialize
cD = [0, sqrt(4*m*k) - 1, sqrt(4*m*k) + 1, sqrt(4*m*k)]; %damping
N = 100;
t = linspace(0,N);
c = 0;
y0 = [1,1];

%Function
f = @(t,y) (-c/m)*y(2) - (k*y(1))/m; %Slope of second dervative
fA = @(t) A*cos(wo*t) + B*sin(wo*t); %Analytic 

%plot
figure(1);clf(1);
plot(t, fA(t)); %Analytic
hold on

%Iterate
for i = 1:4
    c = cD(i);
    f = @(x,u) [u(2)
        -(c/m)*u(2)-(k/m)*u(1)];
    [t,y] = ode45(f,t,y0);
    plot(t,y(:,1));
end

legend("Analytic", "c=0","c=(4*m*k)^(0.5)-1","(4*m*k)^(0.5)+1","(4*m*k)^(0.5)");
xlabel("Time")
ylabel("Position")
%% Problem 2
%Variable Inputs
L = 10;
T = 50;
w = 5;

%initial conditions
u1(1) = 0;
guess = 0;
u2(1) = guess;

%Functions
f = @(x,u) (w/T)*(1+(u2^2))^0.5;

%initialize
N = 100;
x = linspace(0, L, N);
h = x(2) - x(1);
u1 = zeros(1, N); % temp, T
u2 = zeros(1, N); %dT/Dx

%Iterate
    %Loop
    for i = 1:N-1
        %Heun's
        %predictor
        u1_s = u1(i) + h * u2(i);
        u2_s = u2(i) + h * f(x(i), u1(i));
        %corresctor
        u1(i+1) = u1(i) + h/2 * (u2(i) + u2_s);
        u2(i+1) = u2(i) + h/2 * (f(x(i), u1(i)) + f(x(i+1), u1_s));
    end

%plot
figure(2)
plot(x,u1)
xlabel('Distance L')
ylabel('Height')