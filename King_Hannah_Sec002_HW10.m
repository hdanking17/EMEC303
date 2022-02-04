%% EMEC303
% HW 10
% Hannah King

clear; clc;
%% Problem 1
%Number of data points
N = 8;

%Data Values
x = [1, 2, 3, 3.5, 4.25, 5, 7, 10];
y = [0.2, 2.5, 4, 4.1, 4.4, 4, 5, 8.1];

%Plot points
figure(1); clf(1);
plot(x,y,'.');

%Matrix
M = zeros(2,2);
b = zeros(2,1);

M(1,1) = N;
M(1,2) = sum(x);
M(2,1) = sum(x);
M(2,2) = sum(x.^2);

b(1,1) = sum(y);
b(2,1) = sum(y.*x);

a = M\b

ytrend = a(1) + a(2)*x;

hold on
plot(x,ytrend)
hold off

%Sr
Sr = 0;
for i = 1:N
    Sr = Sr + (y(i)-a(1)-a(2)*x(i))^2;
end

%St
St = 0;
for i = 1:N
    St = St + (y(i)-mean(y))^2;
end

%Coeff det
R2 = (St-Sr)/St;
disp("R2: " + R2)

%% Problem 2a
%Number of data points
N = 1000;

%Create N equally spaced data points from 0 to 5
x = linspace(0,5,N);

%Constants
a1 = -3;
a2 = 4;
a3 = 9;

%Data values
y = a1 + a2*x + a3*x.^2 + 5*randn(1,N);

%Plot points
figure(2); clf(2);
plot(x,y,'.');

%Matrix
%Preallocating matrices
f{1} = @(x) 1;
f{2} = @(x) x;
f{3} = @(x) x.^2;

U = length(f);  %Number of basis functions (terms in fit a1+a2*x)
M = zeros(N,U);
b = zeros(N,1);
for i = 1:N
    for j = 1:U
        M(i,j) = f{j}(x(i));
    end
    b(i,1) = y(i);       %b value for the i^th data poit
end

%Solve for a values
a = (M'*M)\(M'*b)

ytrend = a(1) + a(2)*x + a(3)*x.^2;

hold on
plot(x,ytrend)
hold off

%% Problem 2b
% Number of data points
N = 1000;

% Create N equally spaced data points from 0 to 5
x = linspace(0,5,N);

%Constants
a1 = 5;
a2 = -2;
a3 = 5;
a4 = 20;

% Data values
y = a1 + a2*x.^2 + a3*exp(x) + a4*sin(3*pi*x) +...
    5*randn(1,N);

% Plot points
figure(3); clf(3);
plot(x,y,'.');

%Matrix
%Preallocating matrices
f{1} = @(x) 1;
f{2} = @(x) x.^2;
f{3} = @(x) exp(x);
f{4} = @(x) sin(3*pi*x);

U = length(f);  %Number of basis functions (terms in fit a1+a2*x)
M = zeros(N,U);
b = zeros(N,1);
for i = 1:N
    for j = 1:U
        M(i,j) = f{j}(x(i));
    end
    b(i,1) = y(i);       %b value for the i^th data poit
end

%Solve for a values
a =(M'*M)\(M'*b)

ytrend = a(1) + a(2)*x.^2 + a(3)*exp(x) + a(4)*sin(3*pi*x);

hold on
plot(x,ytrend)
hold off