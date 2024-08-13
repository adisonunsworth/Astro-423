%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This program does the necessary calculations for HW04 Astro423
%
%
%
%
%  Author        : C2C Adison Unsworth                   Sept 2024
%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
format long

% Define constants
mu=398600.5;            %Gravitational parameter for Earth

% Calculate the RMOEs and the drift rate of the NMC (𝑦𝑦̇_c). How far does the yc drift in 
% one period? Caution: Note the units given for relative position and velocity

a = 42164; % Semimajor axis km
n = sqrt(mu/a^3); % Mean motion (tgt circular orbit)

x = 2; % km
y = 6; % km
z = 1; % km

xdot0 = 0.0729 /1000; % km/s
ydot0 = -0.2188 /1000; % km/s
zdot0 = 0 /1000; % km/s

% RMOEs

x_c = 4 * x + 2 * ydot0/n;

ydot_c = -1.5 * n * x_c;

y_c = y - 2 * xdot0/n; % km

y_cr = ydot_c * 2 * pi() / n; % km/rev

C = sqrt((x_c - x)^2 + ((y_c - y)/2)^2); % Semi-minor axis (relative)

M = acos((x_c - x) / C);

z_max = sqrt(z^2 + (zdot0 / 2)^2); 

delM = asin(z/z_max);

if xdot0 < 0 
    M = 2*pi - M;
end

if zdot0 < 0
    delM = pi() - delM;
end

if z < 0 && zdot0 > 0 
    delM = 2*pi() + delM; 
end



x_c

y_c

ydot_c

y_cr

C

M = M *180/pi()

z_max 

delM = delM *180/pi()