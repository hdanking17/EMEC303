%% EMEC303
% HW 8
% Hannah King

clear; clc;
%% Problem 1a
%Inputs
L = 0.1; 
D = 1*10^-5;

n=100;

%Initial Conditions
t = 0;
x = linspace(0,L,n);
T = zeros(1,n);
dt = 0.01;
dx= x(2);

for i=1:n
    T(i) = sin((2*pi*x(i))/L);
end

%Iterate
Tnew=T;

for N = 1:1000
    %Time update
    t = t + dt;
    
    for i = 2:n-1
        Tnew(i) = T(i) + dt*(D*(T(i-1) -...
            2*T(i) + T(i+1))/(dx^2));
    end
    T = Tnew;
    
    %Boundary Conditions
    T(1) = 0;
    T(n) = 0;
    
    %Plot
    if rem(N,100)== 0
        figure(1); clf(1)
        hold on
        plot(x,T)
        xlabel("Lenght (x)")
        ylabel("Temp (C)")
        string=sprintf('Time = %7.3f',t);
        title(string)
        hold off
    end
end
