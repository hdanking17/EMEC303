%% EMEC 303 HW6
%  Lance Nichols
%  Section-002
%  9/25/2020

clear; clc;
%% Problem 1: Fluids Value Problem
% * (a) See attached
% * (b) See attached

%% Problem 2: Boundary Value Problem
% * (a) See attached
% * (b) 

% Properties
n = 50;
q = -1000;
L = 10;
E = 200e9;
I = 4e-6;

dx = L/n;

% Instantiate matrices
A = zeros(n);
B = zeros(n,1);

% Add initial conditions
A(1,1) = 1;
A(2,2) = 1;
A(length(A)-1,length(A)-1) = 1;
A(length(A),length(A)) = 1;

for i = 3:n-2
    A(i,i-2) = dx^-4;
    A(i,i+2) = dx^-4;
    A(i,i-1) = -4*dx^-4;
    A(i,i+1) = -4*dx^-4;
    A(i,i) = 6*dx^-4;
    
    B(i,1) = q/(E*I);
end

T=A\B;

plot(linspace(0,L,50),T)
title("Defection of the Beam")
xlabel("Length (m)")
ylabel("deflection (m)")