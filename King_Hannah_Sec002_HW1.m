%% EMEC 303 
%HW1 
%Hannah King

%% Problem 2a
%Create an array
A = zeros(1,20);
A(2) = 1;

%Define Variables
x = 3;

%Loop Over Array
while x<21
    A(x) = A(x-1)+ A(x-2);
    sum = A(x);
    x=x+1;
end

disp('Fibonacci Sequence: ')
disp(sum)

%% Problem 2b
%Plot
figure(1); clf(1)
fplot(@(x) sin(x), [-pi pi], 'b', 1000)
xlabel('x')
ylabel('y')
set(gca, 'Fontsize', 16)
hold on
fplot(@(x) cos(x), [-pi pi], 'r', 1000)
legend('sin(x)', 'cos(x)')

%% Problem 2c
figure(2); clf(2);
imshow('discord.jpg')
title('Roll For It')
xlabel('Life is a game and the dice are weighted against me')