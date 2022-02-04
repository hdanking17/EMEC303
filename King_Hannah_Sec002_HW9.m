%% EMEC303
% HW 9
% Hannah King

clear; clc;
%% Problem 1
N=4; % Number of elements
Dl=2; %Diameter of large cylinder m
Ds=1; %Diameter of small cylinder m
f=2000; % N
L=8; %total length m
E=70*10^9; %Pa or N/m^2 

% Undeformed node locations
x=linspace(0,L,N+1);
dx = x(2) - x(1);

%compute stiffness
kl=pi*(Dl/2)^2*E / dx; %large cylinder spring N/m
ks=pi*(Ds/2)^2*E / dx; %small sylinder spring N/m

%disp(kl)
%disp(ks)

% Prealloate the matrix
K=zeros(N+1,N+1);
F=zeros(N+1,1);

% Fill in matrices
for i=2:N
    K(i,i)=kl + ks;  % Diagonal k(i-1) + k(i)
    K(i,i-1)=-ks;     % Lower diagonal -k(i-1)
    K(i,i+1)=-ks;     % Upper diagonal -k(i)
end
% First row
K(1,1)=1;
% Second row and fourth row
K(2,1)=-kl;
K(N,N+1)=-kl;
% Third row
K(3,3) = ks*2;
F(3,1)=f;
% Last row
K(N+1,N+1)=1;
K(N+1,N)=0;

% Solve for the displacements
U=K\F;

disp(U)
disp(K)

%% Problem 2
N=50; % Number of elements
h=1; %m
L1=0.1; % m
w=1000; %N
E=70*10^9; %Pa or N/m^2
g=9.1; %gravitational acceleration
p= 2710; %kg/m^3 desnity 

% Undeformed node locations
y=linspace(0,L1,N+1);
dy=y(2)-y(1);

min_stress = -1;
TrueL2 = 0;

for L2=y(1):y(i):L1
    d=(L1-L2)/2;
    theta=atan(h/d);
    v=1/3 * ((L1)^2 + L1*L2 + (L2)^2) *h;
    
%stiffness
k=zeros(1,N);
for i=1:N
    dnew=(h-(dy*i))/tan(theta);
    Lx=L2+(2*dnew);
    a=Lx^2;
    k(i)=a*E/dy;
end
% Prealloate the matrix
K=zeros(N+1,N+1);
F=zeros(N+1,1);

% Fill in matrices
for i=2:N
    K(i,i)=k(i-1) + k(i);  % Diagonal k(i-1) + k(i)
    K(i,i-1)=-k(i-1);     % Lower diagonal -k(i-1)
    K(i,i+1)=-k(i);     % Upper diagonal -k(i)
    F(i,1)=v*g*w*p;
end
% First row
K(1,1)=1;
% Last row
K(N+1,N+1)=k(N);
K(N+1,N)=-k(N);
F(N+1,1)=0;

% Solve for the displacements
U=K\F;

% Strain
strain=zeros(1,N);
for i=1:N
    strain(i)=U(i+1)-U(i);
end

% Stress
stress = strain*E;
stress_min = max(stress);

if 0 > min_stress || min_stress > stress_min
    min_stress = stress_min;
    TrueL2 = L2;
end

end

disp(TrueL2)
disp(min_stress)
