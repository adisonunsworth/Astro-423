%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A423_HW_02
%
%  This program calculates initial conditions required to obtain a 
%  periodic relative motion trajectory using the Energy Matching condition 
%  and plots the trajectory obtained by numerically integrating the 
%  General non-linear EOM for satellite ProxOps
%
%  Author        : Lt Col Greg Frey                   Jan 2023
%  Modified by   :
%  Inputs        :
%
%  Outputs       :
%    x0_1 - One initial value for x0 (radial initial position) for periodic motion [km]
%    x0_2 - Second initial value for x0 for periodic motion [km]
%    xout - Numerically integrated state trajectory [various]
%    tout - Time indices for the numerically indicated trajecotry [sec]
%
%  Coupling
%   nonlinEOM.m
%
%  References             :
%    Adapted from Example 4.1 in Spacecraft Formation Flying by Alfriend 
%    et. al., 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%Define Constants
mu=398600.5;    %Gravitational Parameter for Earth

%Define Orbit parameters for the Target spacecraft
a_tgt=8000;             %semi-major axis [km]
e_tgt=0.1;              %eccentricity [-]
n_tgt=sqrt(mu/a_tgt^3); %mean motion [rad/sec]

%% Specify initial conditions (ICs)
%Relative Position and velocity ICs for the chase
y0=0;                               %Initial position (y) [km]
z0=0.1*a_tgt;                       %Initial position (z) [km]
xd0=0.02*a_tgt*n_tgt;               %Initial velocity (xdot) [km/sec]
yd0=0.02*a_tgt*n_tgt;               %Initial velocity (ydot) [km/sec]
zd0=0;                              %Initial velocity (zdot) [km/sec]

%Target ICs 
theta0=0;               %Initial argument of latitude/true anomaly for the target [rad]                   
r0=a_tgt*(1-e_tgt);     %Initial orbital radius for the target (@ perigee) [km]
rd0=0;                  %Initial rate of change of target orbital radius (0 because target is at perigee, a minimum value for r) [km/sec]
thetad0=sqrt(mu/(a_tgt^3*(1-e_tgt^2)^3))*(1+e_tgt*cos(theta0))^2; % Initial rate of change of target arg. of latitude (this is also the angular velocity of the RIC frame wrt the Inertial frame at t=0) [rad/sec]

%% Solve for initial x-position (x0) to yield a 1:1 periodic trajectory
%  using the energy matching condition
% Possible Approach: 
% -Plot the energy matching condition (Ematch) as a function of x0.
% -Use the plot to obtain initial guesses for x0 which yield Emax=0
% -Use a numerical solver to obtain solutions for x0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ****INPUT YOUR CODE HERE****
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x0
eqn =  1/2 *((xd0 - thetad0*y0 + rd0)^2 + (yd0 + thetad0 * (x0 + r0) )^2 + zd0^2) -  (mu)/ (sqrt((r0+x0)^2+y0^2+z0^2)) + mu/(2*a_tgt) == 0; % Equal Energy Equation

% fplot(eqn,x0, [-1000, 1000])

x1 = vpasolve(eqn,x0,[-10000000 0]); % Solve for first equal energy state
x0 = vpasolve(eqn,x0,[-1000 1000]); % Solve for first equal energy state



%% Set Up  Numerical Integration of EOM starting from the ICs obtained above

X=[x0 y0 z0 xd0 yd0 zd0 r0 rd0 theta0 thetad0]'; %Xi is the Initial state vector 
% composed as follows: [x y z xdot ydot zdot r_tgt r_tgtdot theta_tgt theta_tgtdot]'

% Specify time for simulation
tstep=30; %Time step for numerical integration [sec]
P_tgt=2*pi*sqrt(a_tgt^3/mu); %Period of target orbit [seconds]
tf=4*P_tgt; %Final simulation time defined so ~4 orbital periods are simulated [sec]
tspan=0:tstep:tf; %Vector of times for simulation [sec]





% Numerically Integrate using ode45
[tout,xout]=ode45(@nonlinEOM,tspan,double(X)); %% YOU COMPLETE THE FUNCTION nonlinEOM


%% Plot Results
%Create a figure with three plots: x vs. time, y vs. time, z vs. time
figure
subplot(3,1,1)
plot(tout(:,1), xout(:,1))
xlabel('time [sec]')
ylabel('x [km]')
grid on

subplot(3,1,2)
plot(tout(:,1), xout(:,2))
xlabel('time [sec]')
ylabel('y [km]')
grid on

subplot(3,1,3)
plot(tout(:,1), xout(:,3))
xlabel('time [sec]')
ylabel('z [km]')
grid on

%Create a plot containing a 3-dimensional view of the position trajectory
figure
plot3(xout(:,1),xout(:,2),xout(:,3))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on

