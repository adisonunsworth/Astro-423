function [dV1, dV2, Xf]=HCW_Targeting(r_tgt, X0, Xfd, TOF)
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
%    r_tgt - target orbit radius [km]
%    X0 - Initial RIC-frame state vector for chase [x0 y0 z0 xdot0 ydot0 zdot0]', units of km and km/s
%    Xfd - Desired final RIC-frame state vector for chase [xf yf zf xdotf ydotf zdotf]', units of km and km/s
%    TOF - Time of flight to get to desired final state [sec]
%  Outputs       :
%    dV1 - First delta V [km/s]
%    dV2 - Second deltaV [km/s]
%    Xf  - Final chase position and velocity after maneuvers [km, km/s]
%  Coupling
%   HCW.m - Propagates initial state forward in time using the HCW solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Target vehicle information
h_tgt=r_tgt-6378.137;         %Target altitude [km]
n=sqrt(398600.5/r_tgt^3);     %Target mean motion [rad/sec]


%Extract Initial position from the initial state vector [km]
x0=X0(1);
y0=X0(2);
z0=X0(3);
r0=[x0; y0; z0];

%Extract Initial velocity from the initial state vector [km/s]
xd0=X0(4);
yd0=X0(5);
zd0=X0(6);
v0=[xd0; yd0; zd0];

%Extract Desired final position from final state vector
xf=Xfd(1);
yf=Xfd(2);
zf=Xfd(3);
rf=[xf;yf;zf];

%Desired final velocity
xdf=Xfd(4);
ydf=Xfd(5);
zdf=Xfd(6);
vf=[xdf; ydf; zdf];


%% Calculate the two delta Vs required using HCW targeting
%Define sin and cosine terms to be used multiple times in calculations
c=cos(n*TOF);
s=sin(n*TOF);

%Calculate the Phirr and Phivr matrices
Phirr=  [4-3*c      0   0; 
        6*s-6*n*TOF   1   0; 
        0           0   c];
Phivr=  [s/n        (2/n)*(1-c)     0;
        (2/n)*(c-1) -3*TOF+(4*s/n)    0;
        0           0               s/n];

%Calculate the velocity needed at t=0 to reach the desired final position
%at the desired final time. 
vneed=inv(Phivr)*(rf-Phirr*r0);

%Calculate delta V1
dV1=vneed-v0;

%Calculate the position and velocity after the Time of flight
Xf=HCW(h_tgt, [r0;vneed],TOF);  %% ** Use your HCW Function!

%Calculate second dv to obtain the desired final velocity
dV2=vf-Xf(4:6);

%Final state after maneuver
Xf=Xf+[0;0;0;dV2];
end



