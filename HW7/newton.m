function E = newton(mean,ecc)
%===========================================================
%
% Use: [E] = newtonr(mean,ecc)
%
% The function determines the Eccentric anomaly from the mean anomaly and
% the eccentricity. Iterates  until 10^-9 is the divergence
%
% Author: Adison Everett Unsworth  USAFA/20 803-448-3527
%
% Inputs:
%   mean    -   Mean Anomaly            -   Radians
%   ecc     -   Eccentricity            -   Unitless
%
% Outputs:
%   E       -   Eccentric anomaly       -   Radians
%
% Locals:
%
% Globals: None
%
% Constants: None
%
% Coupling: None
%
% Assumptions: None
% 
% References:
%   Astro 310 Equation sheet 2023
%   Fundementals of Astrodynamics and applications fourth edition Vallado page 65
%
% Validated with ASTRO 310 LSN 15 Predicting Orbits Youtube
%   4.367759837929274 Radians
%   mean    = 4.55596
%   ecc     = 0.2
%===========================================================
% create placeholder variable
Eold = mean;
if pi < mean
    Enew = mean - ecc;
else
    Enew = mean + ecc;
end

while abs(Enew - Eold) > 10^-9
    Eold = Enew;
    Enew = Eold + (mean - Eold + ecc*sin(Eold))/(1- ecc*cos(Eold));   
end

E = Enew;

