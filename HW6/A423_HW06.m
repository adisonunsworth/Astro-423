%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HCW_Targeting_HW06
%
%  This program calculates the two delta-Vs required to maneuver a chase 
%  satellite from a specified initial position/velocity to a specified 
%  final position/velocity in a specified amount of time using the HCW 
%  targeting technique.
%
%  Author        : Lt Col Greg Frey                   Feb 2023
%  Modified by   :
%  Inputs        :
%
%  Outputs       :
%    dV1 - First delta V [km/s]
%    dV2 - Second deltaV [km/s]
%    Xf  - Final chase position and velocity after maneuvers [km, km/s]
%  Coupling
%  HCW_Targeting.m - Calculates delta-Vs using HCW Targeting technique   
%  HCW.m - Propagates initial state forward in time using HCW solution
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long

%Target vehicle information
r_tgt=42164; %Geo orbit radius [km]

%%**************YOU COMPLETE THE FOLLOWING CODE****************************
% Fill in the following ICs for the Chase %%
%Initial position [km]
x0= -14.142;
y0= 28.284;
z0= 10;
r0=[x0; y0; z0];

%Initial velocity [km/s]
xd0= .0010313;
yd0= .0020625;
zd0= 0;
v0=[xd0; yd0; zd0];

%Initial state vector
X0=[r0; v0];

% Fill in the following information about the desired trajectory %%
TOF= 3600*11; %Time from now until reaching desired location [seconds] 

%Desired final position
xf= 0;
yf= 25;
zf= 0;
rf=[xf;yf;zf];

%Desired final velocity
xdf= 0;
ydf= 0;
zdf= 0;
vf=[xdf; ydf; zdf];

%Desired final state vector
Xfd=[rf; vf];
%%*************************************************************************

%% Calculate the two delta Vs required using HCW targeting
[dV1, dV2, Xf]=HCW_Targeting(r_tgt, X0, Xfd, TOF); %Note that HCW_Targeting called HCW.m. Use your HCW.m from HW03!

%Output results
fprintf('Delta V_1 in the RIC frame [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n', dV1(1), dV1(2), dV1(3))
fprintf('Delta V_2 in the RIC frame [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n', dV2(1), dV2(2), dV2(3))
fprintf('Final Position in the RIC frame [km]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n', Xf(1), Xf(2), Xf(3))
fprintf('Final Velocity in the RIC frame [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n', Xf(4), Xf(5), Xf(6))
