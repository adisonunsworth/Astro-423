function wgs84=wgs84data
%==========================================================================
%
% use: wgs84 = wgs84data
%                
% This function creates a structure that contains wgs84 values and several
% other constants.
%
% Author
%  Originally written by Capt Dave Vallado
%  Modified and Extended for Ada by Dr Ron Lisowski
%  Extended from DFASMath.adb by Thomas L. Yoder, LtCol, Spring 00
%  Made WGS84 and derived values consistent by Scott R. Dahlke, Spring 08
%  Modified by LtCol Luke Sauter into a structured array 2017
%  Added Documentation Scott Dahlke, 2023-06-20
%
% Inputs
%   None
%
% Outputs
%   wgs84 - structure with wgs84 values and other constants 
%           wgs84.Deg         radians to degrees factor deg/rad
%           wgs84.Rad         degrees to radians factor rad/deg
%           wgs84.MU          *** gravitational parameter of Earth km^3/sec^2
%           wgs84.RE          *** radius of Earth km
%           wgs84.OmegaEarth  *** rotation rate of Earth rad/sec
%           wgs84.SidePerSol  Sidereal Days/Solar Day
%           wgs84.RadPerDay   rotation rate of Earth rad/day
%           wgs84.SecDay      sec/day
%           wgs84.Flat        *** flattening of the earth (unitless)
%           wgs84.EEsqrd      eccentricity of Earth's ellipsoid squared (unitless)
%           wgs84.EEarth      eccentricity of Earth's ellipsoid (unitless)
%           wgs84.J2          J2 geopotential of Earth (oblateness) (unitless)
%           wgs84.J3          J3 geopotential of Earth (unitless)
%           wgs84.J4          J4 geopotential of Earth (unitless) 
%           wgs84.GMM         gravitational parameter of the Moon km^3/sec^2
%           wgs84.GMS         gravitational parameter of the Sun  km^3/sec^2
%           wgs84.AU          distance from Sun to Earth km
%           wgs84.HalfPI      pi/2 unitless
%           wgs84.TwoPI       2*pi unitless
%           wgs84.Zero_IE     Small number for incl & ecc purposes
%           wgs84.Small       Small number used for tolerance purposes
%           *** indicates the four acutal WGS84 values 
%
% References
%   https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84
%
%==========================================================================
   
% Build Structured array
wgs84.Deg        = 180.0/pi;                   % radians to degrees factor deg/rad
wgs84.Rad        = pi/180.0;                   % degrees to radians factor rad/deg
wgs84.MU         = 398600.5;                   % *** gravitational parameter of Earth km^3/sec^2
wgs84.RE         = 6378.137;                   % *** radius of Earth km
wgs84.OmegaEarth = 0.00007292115;              % *** rotation rate of Earth rad/sec
wgs84.SidePerSol = 0.00007292115*86400/(2*pi); % Sidereal Days/Solar Day
wgs84.RadPerDay  = 0.00007292115*86400;        % rotation rate of Earth rad/day
wgs84.SecDay     = 86400.0;                    % sec/day
wgs84.Flat       = 1.0/298.257223563;          % *** flattening of the earth (unitless)
wgs84.EEsqrd     = (2.0-(1.0/298.257223563))*(1.0/298.257223563);       % eccentricity of Earth's ellipsoid squared (unitless)
wgs84.EEarth     = sqrt((2.0-(1.0/298.257223563))*(1.0/298.257223563)); % eccentricity of Earth's ellipsoid (unitless)
wgs84.J2         = 0.00108263;                 % J2 geopotential of Earth (oblateness) (unitless)
wgs84.J3         = -0.00000254;                % J3 geopotential of Earth (unitless)
wgs84.J4         = -0.00000161;                % J4 geopotential of Earth (unitless) 
wgs84.GMM        = 4902.774191985;             % gravitational parameter of the Moon km^3/sec^2
wgs84.GMS        = 1.32712438E11;              % gravitational parameter of the Sun  km^3/sec^2
wgs84.AU         = 149597870.0;                % distance from Sun to Earth km
wgs84.HalfPI     = pi/2.0;                     % pi/2 unitless
wgs84.TwoPI      = 2.0*pi;                     % 2*pi unitless
wgs84.Zero_IE    = 0.015;                      % Small number for incl & ecc purposes
wgs84.Small      = 1.0E-6;                     % Small number used for tolerance purposes

end         
