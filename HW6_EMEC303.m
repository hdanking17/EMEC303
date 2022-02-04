%% EMEC303
% HW 6
% Hannah King

clear; clc;
%% Problem 2
%inputs
L = 10; %m
q = -1000; %N/m
E = 200*10^9; %Pa
I = 4*10^-6; %m^4

%Matrix method
% initialize matrices
N=50; % number of discrete points on object (meshsize in 2D)
A = zeros(N,N);
b = zeros(N,1);
dx = L/N; %step size
xm = linspace(0,L,N);
RHS = q/(E*I);

% BC
A(1,1) = 1;
A(2,2) = 1;
A(N-1, N-1) = 1;
A(N,N) = 1;

% Fill matrices
for i=3:N-2
    A(i,i) = -4;
    A(i,i-2) = 1;
    A(i,i-1) = 1;
    A(i,i+1) = 1;
    A(i,i+2) = 1;
    b(i) = RHS;
end

% Solve
T_mat = -A\b;

%plot
figure(1); clf(1);
plot(xm,T_mat) % Matrix method
title('Deflection of the Beam')
xlabel('Length')
ylabel('Deflection')