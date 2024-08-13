function Xf = HCW(StartAlt,Xi,t)
% Given an altitude of the target, initial conditions of the interceptor
% and the time, this function calculates the target's state

% INPUTS:  StartAlt [km]
%          Xi [km; km; km; km/s;  km/s;  km/s ] in RIC frame
%             [x0; y0; z0; xdot0; ydot0; zdot0]
%          t [s]
%
% OUTPUTS: Xf [km; km; km; km/s; km/s; km/s] in RIC frame
%             [x;  y;  z;  xdot; ydot; zdot]
% define global constants
%global MU RE
MU=398600.5; 						%% km^3/sec^2
RE=6378.137; 						%% km
%% write HCW equations here...
x0 = Xi(1);
y0 = Xi(2);
z0 = Xi(3);
xdot0 = Xi(4);
ydot0 = Xi(5);
zdot0 = Xi(6);

n = sqrt(MU/((RE+StartAlt)^3));

x = (xdot0/n) *sin(n*t) - (3*x0 +2*ydot0/n) * cos(n*t) +(4*x0 + 2 * ydot0/n);
y = (6*x0 + 4*ydot0/n) *sin(n*t) + 2 * xdot0/n *cos(n*t) - (6*n*x0 + 3*ydot0)*t + (y0 - 2*xdot0/n);
z = z0 * cos(n*t) + zdot0/n * sin(n*t);
xdot = xdot0 * cos(n*t) + (3 * n*x0 + 2*ydot0) *sin(n*t);
ydot = (6*n*x0 + 4*ydot0) * cos(n*t) - 2*xdot0 * sin(n*t) - (6*n*x0+3*ydot0);
zdot = -z0 * n * sin(n * t) + zdot0 *cos(n*t);
Xf = [x;  y;  z;  xdot; ydot; zdot];

end