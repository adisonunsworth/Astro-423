function [rpqw,vpqw]=RVpqw(a,ecc,nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use: [rpqw,vpqw]=RVpqw(a,ecc,nu)
%
% The function "RVpqw" computes R and V in the pqw frame
%
% Author: Scott Dahlke  USAFA/DFAS  719-333-4110     1 Jan 09
%
% Inputs:
%   a - semimajor axis (km)
%   ecc - eccentricity
%   nu - true anomaly (rad) 
%
% Outputs:
%   rpqw - position in the pqw frame (km)
%   vpqw - velocity in the pqw frame (km/s)
%
% Globals: MU
%
% Constants: None
%
% Coupling: None
% 
% References:
%   COE's to RV Lesson of Astro 201
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global MU

% compute the parameter p
p = a*(1-ecc^2);

% compute the position vector
rpqw = p/(1 + ecc*cos(nu))*[cos(nu);sin(nu);0];

% compute the velocity vector
vpqw = sqrt(MU/p).*[-sin(nu);(ecc+cos(nu));0];
    