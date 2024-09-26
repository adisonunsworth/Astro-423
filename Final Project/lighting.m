function [CATS,local_time] = lighting(zulu_time,longitude,rho_tgt_chase,DOY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [CATS,local_time] = lighting(zulu_time,longitude,rho_tgt_chase,DOY)
%
%
% Inputs:
%           zulu_time:      hr
%           longitude:      degrees
%           rho_tgt_chase:  3x1 vector [x,y,z] 
%           rho_tgt_chase:  3x1 vector [x,y,z] km x is virtical, y is
%                           horizontal and negative and z is normal
%           DOY:            
% Outputs:
%           CATS:           Angle between sun and tgt and chase degrees, 
%                           closer to 0 is better.
%           local_time:     Time at location    - hours
%
%
% Author:   C2C Adison Unsworth 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes DOY based on Changes in zulu time 
while zulu_time >= 24 
    if zulu_time >= 24 
        DOY = DOY +1;
        zulu_time = zulu_time - 24;
    end
end

zulu_offset = longitude / 15; % hrs

local_time  = zulu_time + zulu_offset; % hrs


t = local_time * 3600; % time since local midnight - sec

t = local_time * 86400; % time since local midnight - sec


theta = 360 * t / (86400); % angle the sun has traveled since local midnight

beta = -23.45 * cosd(360 * (DOY + 10)/(365)); % angle from r_sun_hat to tgt orbital plane

r_sun_xy = 1/sqrt((1 + (tand(beta))^2));

r_sun_hat = [-r_sun_xy * cosd(theta); r_sun_xy * sind(theta)  ;r_sun_xy * tand(beta)];

CATS  = acosd(dot(r_sun_hat,rho_tgt_chase)/(norm(rho_tgt_chase)));
