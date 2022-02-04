%% EMEC303
% HW 7
% Hannah King

clear; clc;
%% Problem 1b
% Domain size
Lx=1;
Ly=1;

% Temps on the Boundaries
Tr=10;
Tl=10;
Tb=10;
Tt=10;

% Number of grid points
Nx=20;
Ny=20;

% Create grid
x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);

% Build matrics
A=zeros(Nx*Ny,Nx*Ny);
B=zeros(Nx*Ny,1);

% Fill in rows for interior grid points
for i=2:Nx-1
    for j=2:Ny-1
        % Compute the row number n
        n=i+(j-1)*Nx;
        % Put in non-zero terms
        A(n,n)=-4;
        A(n,n-1)= 1;
        A(n,n+1)= 1;
        A(n,n-Nx)= 1;
        A(n,n+Nx)= 1;
        Actr=1;
        sig=0.1;
        B(n,1)=  - (Actr*exp(-((x(i)-.5)^2 + (y(j)-.5)^2) / (2*sig^2)));
    end
end

% Left boundary
i=1;
for j=1:Ny
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=10;
end

% Right boundary
i=Nx;
for j=1:Ny
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=10;
end

% Bottom boundary
j=1;
for i=2:Nx-1 % Skipping corners
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=10;
end

% Top boundary
j=Ny;
for i=2:Nx-1 % Skipping corners
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=10;
end

% Solve for the temperature
Tv=A\B;  

% Convert 1D vector into 2D array of temps
h=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        % Compute the row number n
        n=i+(j-1)*Nx;
        % Put Tv(n) into T(i,j)
        h(i,j)=Tv(n);
    end
end

% Plot
figure(1); clf(1)
surf(x,y,h')
xlabel('x (m)')
ylabel('y (m)')
zlabel('T(x,y) in Celcius')
title('Temperature of Plate Heated From Center Cooled on all Edges')

%% Problem 1c
% Domain size
Lx=1;
Ly=1;

% Temps on the Boundaries
Tr=10;
Tl=10;
Tb=10;
Tt=10;

% Number of grid points
Nx=20;
Ny=20;

% Create grid
x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);

% Build matrics
A=zeros(Nx*Ny,Nx*Ny);
B=zeros(Nx*Ny,1);

% Fill in rows for interior grid points
for i=2:Nx-1
    for j=2:Ny-1
        % Compute the row number n
        n=i+(j-1)*Nx;
        % Put in non-zero terms
        A(n,n)=-4;
        A(n,n-1)= 1;
        A(n,n+1)= 1;
        A(n,n-Nx)= 1;
        A(n,n+Nx)= 1;
        Actr=1;
        sig=0.1;
        B(n,1)=  - (Actr*exp(-((x(i)-.5)^2 + (y(j)-.5)^2) / (2*sig^2)));
    end
end

% Left boundary
i=1;
for j=1:Ny
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=0;
end

% Right boundary
i=Nx;
for j=1:Ny
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=0;
end

% Bottom boundary
j=1;
for i=2:Nx-1 % Skipping corners
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=10;
end

% Top boundary
j=Ny;
for i=2:Nx-1 % Skipping corners
    % Compute the row number n
    n=i+(j-1)*Nx;
    % Put in non-zero terms
    A(n,n)=1;
    B(n,1)=10;
end

% Solve for the temperature
Tv=A\B;  

% Convert 1D vector into 2D array of temps
h=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        % Compute the row number n
        n=i+(j-1)*Nx;
        % Put Tv(n) into T(i,j)
        h(i,j)=Tv(n);
    end
end

% Plot
figure(2); clf(2)
surf(x,y,h')
xlabel('x (m)')
ylabel('y (m)')
zlabel('T(x,y) in Celcius')
title('Temperature of Plate Heated From Center Cooled on Top and Bottom Edges')

%% Problem 1d
%imputs
Lx = 1; %size of domain
Ly = 1;
Nx = 20; %Number of grid points
Ny = 20;
sigma = 0.1;
Actr = 1;

Niter = 10000; %max number of iterations to allow the results to converge
tol = 1e-6; %convergence tolerance

%grid
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
dx = x(2) - x(1);
dy = y(2)-y(1);

%initial guess for solution
Phi = zeros(Nx, Ny);

%iterate
PhiNew = Phi; %preallocate

for n = 1:Niter
    %loop over interior points
    for i = 2:Nx-1
        for j = 2:Ny-1
            PhiNew(i,j) = (Actr*exp(-((x(i)-.5)^2+(y(j)-.5)^2)/...
                (2*sigma^2))) + (0.25 * (Phi(i-1,j)+Phi(i+1,j) ...
                + Phi(i,j-1) + Phi(i,j+1)));
        end
    end
    
    %Boundary Conditions
    %left
    i = 1;
    for j = 1:Ny
        PhiNew(i,j) = 10;
    end

    %right
    i = Nx;
    for j = 1:Ny
        PhiNew(i,j) = 10;
    end

    %bottom
    j = 1;
    for i = 2:Nx-1
        PhiNew(i,j) = 10;
    end

    %top
    j = Ny;
    for i = 2:Ny-1
        PhiNew(i,j) = 10;
    end

%check if converged
residual = max(abs(PhiNew(:) - Phi(:)));
if residual < tol
    break
end

Phi = PhiNew;

%plot
if rem(n,100) == 0
    figure(3);clf(3)
   surf(x,y,Phi')
   xlabel('x (m)')
   ylabel('y (m)')
   zlabel('Phi(x,y) in Celcius')
   title('Temperature of Plate Heated From Center - Liebmanns')
   set(gca,'Fontsize', 20)
   drawnow
end
end