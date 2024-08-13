function [rijk,vijk]=RVijk(rpqw,vpqw,incl,raan,argp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use: [rijk,vijk]=RVijk(rpqw,vpqw,incl,raan,argp)
%
% The function "RVijk" rotates R and V in the pqw frame to 
% the ijk frame
%
% Author: Scott Dahlke  USAFA/DFAS  719-333-4110    1 Jan 09
%
% Inputs:
%   rpqw - position in the pqw frame (km)
%   vpqw - velocity in the pqw frame (km/s)
%   incl - inclination (rad)
%   raan - right ascension of the ascending node (rad)
%   argp - argument of perigee (rad)
%
% Outputs:
%   rijk - position in the ijk frame (km)
%   vijk - velocity in the ijk frame (km/s)
%
% Globals: None
%
% Constants: None
%
% Coupling: axisrot
% 
% References:
%   COE's to RV Lesson of Astro 201
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rotate the position vector
rtemp1 = axisrot(rpqw,3,-argp);
rtemp2 = axisrot(rtemp1,1,-incl);
rijk   = axisrot(rtemp2,3,-raan);

% rotate the velocity vector
vtemp1 = axisrot(vpqw,3,-argp);
vtemp2 = axisrot(vtemp1,1,-incl);
vijk   = axisrot(vtemp2,3,-raan);

