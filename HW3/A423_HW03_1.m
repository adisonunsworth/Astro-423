%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A423_HW_3
%
%  This program calculates initial relative position and velocity based on
%  a scenario where target and chase satellites both have the same semi-major 
%  axis, the target's orbit is circular, the chase's orbit is elliptical 
%  and the chase begins a distance of x km above the target. Using these 
%  initial conditions, the program:
%      - Calculates and plots the value of the Energy Matching Condition 
%        and the HCW condition for periodic motion as a function of x
%      - Calculates and plots the relative motion trajectory for x=5km and 
%        x=500km using the general non-linear EOM for ProxOps and the 
%        HCW solution
%
%  Author        : Lt Col Greg Frey                   Aug 2023
%  Modified by   : C2C Adison Unsworth
%
%  Coupling
%   nonlinEOM.m - General non-linear EOM for proxops
%   ode45 - Matlab's built-in numerical integrator
%   HCW.m - Calculates relative position and velocity using HCW solution
%
%  References             : None
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long

% Define constants
mu=398600.5;            %Gravitational parameter for Earth

% Define target orbital parameters
r_tgt=7000;              %Target orbital radius [km]
n_tgt=sqrt(mu/r_tgt^3);   %Target mean motion [rad/sec]

%% HW03 Part c: 
% Calculate the HCW condition for Periodic Motion and the Energy Matching 
% condition for specified initial conditions

% Define range of x values to be used.
xrange=0:5:500;  %Given range of initial x position [km]

% Initialize vectors to store the values for the HCW condition for periodic
% motion and the Energy matching condition for each x value.
HCW_condition=inf(1, length(xrange));
Ematch=inf(1, length(xrange));

% In loop below: For each x value, calculate the other initial position and 
% velocity components (y, z, xdot, ydot, zdot) as well as the, HCW condition
% and Energy Matching condition. Store values.
counter=1; %Initialize counter for proper indexing.
for x=xrange
    
    %%%%%%%%%%%%%%%%YOU COMPLETE THE FOLLOWING CODE%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Code expressions for initial relative position and initial
    % relative velocity (with respect to the RIC frame). Note that some of
    % these components will depend on the value of x.
    v_tgt   = sqrt(mu/r_tgt);

    v_chase = sqrt(2 *(mu/(r_tgt+x) - mu/(2*r_tgt)));
    
    ydot_I =  v_chase - v_tgt;

    n = sqrt(mu/r_tgt^3);
    
    y = 0 ; %initial in-track position
    z = 0 ; %initial cross-track position
    
    xdot_ijk = 0; %initial radial velocity (wrt the IJK resolved in the RIC frame)
    ydot_ijk = ydot_I; %initial in-track velocity (wrt the IJK resolved in the RIC frame)
    zdot_ijk = 0; %initial cross-track velocity (wrt the IJK resolved in the RIC frame)
    
    

    V = [xdot_ijk; ydot_ijk; zdot_ijk] + cross([0;0;-n], [x;0;0]);
    
    xdot = V(1);
    ydot = V(2);
    zdot = V(3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Using the initial relative position and velocity above, calculate the
    %value of the HCW condition for periodic motion and the Energy Matching
    %condition
    HCWcase=ydot+2*n_tgt*x;   %HCW condition for periodic motion
    Ematchcase=0.5*((xdot-n_tgt*y)^2+(ydot+n_tgt*(x+r_tgt))^2+zdot^2)-mu*((r_tgt+x)^2+y^2+z^2).^(-0.5)+(mu/(2*r_tgt)); %Energy matching condition
    
    %Store these values
    HCW_condition(counter)=HCWcase;
    Ematch(counter)=Ematchcase;
    
    %Increment counter
    counter=counter+1;
end

%Plot the HCW Condition and the Energy Matching Condition as functions of x
plot(xrange,HCW_condition)
hold on
grid on
plot(xrange,Ematch)
xlabel('x [km]')
legend('HCW Condition','Energy Matching Condition')

%% HW03 Part d: 
% For x=5km and x=500km, simulate the relative motion for one period of the
% target's orbit using the general non-linear EOM for ProxOps and the HCW
% solution.

% Set up time for simulations
P_tgt=2*pi*sqrt(r_tgt^3/mu); %Period of target orbit [sec]
tf=P_tgt; %Final time defined so ~1 orbital period is simulated [sec]
tstep=30; %time step of 30 seconds
tspan=0:tstep:tf; %Vector of times for simulation

for x=[5 500]  %Simulate Relative Motion with x=5km and x=500km
    
    %%%%%%%%%%%%%%%%YOU COMPLETE THE FOLLOWING CODE%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the initial relative position and initial
    % relative velocity (with respect to the RIC frame). Note that some of
    % these components will depend on the value of x.
        
      

    v_tgt   = sqrt(mu/r_tgt);

    v_chase = sqrt(2 *(mu/(r_tgt+x) - mu/(2*r_tgt)));
    
    ydot_I =  v_chase - v_tgt;

    n = sqrt(mu/r_tgt^3);
    
    y = 0 ; %initial in-track position
    z = 0 ; %initial cross-track position
    
    xdot_ijk = 0; %initial radial velocity (wrt the IJK resolved in the RIC frame)
    ydot_ijk = ydot_I; %initial in-track velocity (wrt the IJK resolved in the RIC frame)
    zdot_ijk = 0; %initial cross-track velocity (wrt the IJK resolved in the RIC frame)
    


    V = [xdot_ijk; ydot_ijk; zdot_ijk] + cross([0;0;-n], [x;0;0]);
    
    xdot = V(1);
    ydot = V(2);
    zdot = V(3);
   
    % EOM
    r_tgtdot = 0; %rate of change of target orbital radius [km/s]
    theta_tgt = 0; %initial argument of latitude for target [rad]
    theta_tgtdot = sqrt(mu/r_tgt^3); %initial rate of change of tgt arg of latitude [rad/s]
    
    %Form the intial state vector for non-linear EOM
    Xi_numint=[x y z xdot ydot zdot r_tgt r_tgtdot theta_tgt theta_tgtdot]'; %Initial condition for numerical integration.
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Generate trajectory by numerically integrating the non-linear EOM
    %Integrate general nonlinear EOM using ode45
    
    
    [tout,xout]=ode45(@nonlinEOM,tspan,Xi_numint); %%%% USE YOUR nonlinEOM FUNCTION FROM HW02

    % Generate trajectory using the solution to the HCW equations
    % Define initial state vector for HCW calculations
    Xi_HCW=[x y z xdot ydot zdot]';   %Initial condition for HCW
    
    counter=1;                  %Initialize a counter for proper indexing
    X_HCW=inf(length(tspan),6); %Pre-allocate matrix to store HCW trajectory
    for ii=1:length(tspan)
        X_HCW(counter,:) = HCW(r_tgt-6378.137,Xi_HCW,tspan(ii));  %%%% YOU WRITE THE FUNCTION HCW.m
        counter = counter+1;
    end

    % Plot the relative motion trajectories in the x-y plane
    figure
    plot(xout(:,2),xout(:,1))
    hold on
    plot(X_HCW(:,2),X_HCW(:,1),'--')
    xlabel('y - in-track [km]')
    ylabel('x - radial [km]')
    legend('Numerically Integrated Trajectory (ode45)','HCW Trajectory','Location','southwest')
    if x==5
        title('Relative Trajectory in x-y plane: x=5 km')
    else
        title('Relative Trajectory in x-y plane: x=500 km')
    end
    grid on

end
