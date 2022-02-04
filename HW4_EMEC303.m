%% EMEC303
%HW 4
% Hannah King

clear; clc;
%% Problem 1
% Define inputs
Dt = 1; %m
%D0 = 0.01; %m
D0 = 0.1; %m
g = 9.8; %m/s^2
h = 0.5; %step size
i = 1;

%Define function
f = @(t,y) -D0^2 / Dt^2 * sqrt(2*g*y);

%initial conditions
t(1) = (0);
y(1) = 1; %m

%iterate
while y >= 0
    t(i+1) = t(i) + h;
    ys = y(i) + h*f(t(i), y(i)); %predictor
    y(i+1) = y(i) + h/2 * (f(t(i), y(i)) + f(t(i+1), ys)); %corrector
    i = i+1;
end

figure(1);
clf(1);
plot(t, y, '-')
xlabel('Time (s)')
ylabel('Height (m)')

%% Problem 2
%define inputs
xL = 0;
%xU = 1;
xU = 0.02;
%h = 0.0015; %step size
h = 0.0001; %step size
y0 = 0;

%define function
f = @(x,y) -1000*y + 3000 - 2000*exp(-x);

%arrays
x = [xL:h:xU];
y = zeros(1, length(x));
ye = y;

%initial conditions
y(1) = y0;
ye(1) = y0;

%analytical equation
fA = @(x) 3 - 0.998*exp(-1000*x) - 2.002*exp(-x);

%iterate
for i=1:length(x)-1
    %Euler's
    ye(i+1) = ye(i) + h*f(x(i),ye(i));
    %Runge-Kutta
    k1 = f(x(i), y(i));
    k2 = f(x(i) + 1/2 * h, y(i) + 1/2 * h * k1);
    k3 = f(x(i) + 1/2 * h, y(i) + 1/2 * h * k2);
    k4 = f(x(i) + h, y(i) + k3 * h);
    y(i+1) = y(i) + 1/6 * (k1 + 2*k2 + 2*k3 + k4)*h;
    
    
end

%plot
figure(2);
clf(2);
plot(x, y, 'b-');
hold on;
plot(x, ye, 'k*', 'markersize', 2);
plot(x, fA(x), 'ro', 'markersize', 2);
legend('Runge-Kutta Method', 'Euler''s method', 'Analytic')