%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Phase_1_Master
%
% 1. Your mission planning starts when the X-37C arrives in Geostationary orbit at 1500Z on 21 March 2025. 
%   The X-37C arrives 50 km ahead of the refueling hub, which is stationed over 75o
%   West Longitude. You may use the HCW model for your mission planning. You may assume one day and one period of a 
%   Geostationary orbit are exactly 24 hrs long for your calculations. 
% 
% 2. Calculate ΔVs to enter a 1km x 2km NMC (C=0.5 km) for initial inspection using visual sensors, centered 
%   at the refueling hub. 
%       a. The NMC should have Zmax = 100 m, with the position of Zmax in the NMC chosen to maximize the 
%           safety of the initial NMC.
%       b. The maneuvers to accomplish the ingress and NMC should be designed to: 
%           i. Allow the X-37C to enter the NMC no later than 2.5 days after arriving in GEO
%           ii. Ensure appropriate lighting for inspection with a visual sensor (CATS near 0o within the 
%               X-Y plane). 
%           iii. Utilize less total ΔV if multiple maneuver options are feasible*
%       c. The X-37C must remain outside of a spherical keep out zone centered on the refueling hub with 
%           a radius of 400 m at all times during this initial phase of the ingress and inspection. 
% 
% 3. The X-37C must complete two complete revolutions in the initial inspection NMC to ensure sufficient 
%   time for inspection.
% 
% 4. The initial inspection reveals that the docking port on the +y face of the refueling hub may have been 
%   damaged by the micrometeoroid strike. Calculate delta-Vs needed to modify your trajectory to a 
%   0.1x0.2km NMC (C=0.05 km) centered on the refueling hub with Zmax=0 km to obtain a closer look at the 
%   +y face using an on-board LIDAR system. 
%       a. During this phase of the mission, the X-37C must remain outside of a spherical keep out zone 
%           centered in the refueling hub with a radius of 40 m at all times. 
% 
% 5. After two complete revolutions in the new NMC, data from the secondary inspection reveals that all is 
%   well. Calculate ΔVs needed to ingress and park at 50m in front of the refueling hub, at which point an 
%   automated system will take over to complete the refueling and docking. 
%       a. During this phase of the mission, the X-37C must remain outside of a spherical keep out zone 
%           centered in the refueling hub with a radius of 40 m at all times. 
% 
% 6. Ensure you report the total amount of ΔV used in the ProxOps portion of the mission. 
% 
% * Consider if multiple times of flight will yield the desired trajectory given the constraints and calculate a few 
% options to try to reduce delta-V. You do not have to perform any type of optimization or in-depth analysis to 
% minimize delta-V.
%
%  Author        : C2C Adison Unsworth                   Sept 2024
%  
%
%  Coupling
%   wgs84
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

%Constants
wgs84 = wgs84data;

MU = wgs84.MU; % Gravitational Constant - km^3/sec^2

% Geostationary_Orbit
LOD = 24*3600; % length of day - seconds

period = 24*3600; % Period - seconds -  in the prompt geo orbit is defined by its period
a = ((period/(2*pi()))^2*MU)^(1/3); % Semi-Major Axis - km
r_tgt = a; % due to the circular nature of the orbit
n = sqrt(MU/a^3); %mean motion

% Other Inputs
    
    zulu_time_i = 15; % Hour
    day_of_year = 80; %21 March 2025
    y_ci = 50; % Chase begins 50 km ahead of tgt - km
    y_cf = 0; % final desired center of NMC
    z_max0 = 0; % inital Condition
    longitude_i = -75; % initial longitude degrees west of zulu;
    

    % Desired Characteristics
    z_max1 = .1; % maximum distance in the z direction km
    C = 0.5; % semi-minor axis km 
    x_c0 = 0; % x center of the ellipse km

    rho_tgt_chase_i = [0;y_ci;0]; % inital position of the tgt to the chase - km
    rho_dot_tgt_chase_i = [0;0;0]; % inital velocity of the tgt to the chase - km/s
    
    rho_tgt_chase_f = [0;C*2;0]; % position at final burn
    
% Inital Cats
[CATS_i,local_time_i] = lighting(zulu_time_i,longitude_i,rho_tgt_chase_i,day_of_year);

% Desired Transfer Characteristics
max_time = 2.5 * LOD; % max time to get into an orbit defined by C - seconds
keep_out = .600; % distance to stay away from tgt -  km

% determine Delta v needed to enter into NMC from stationary y = C*2 km/s
y_c1 = C*2; % km/s

delv_s_to_NMC = [-n/2 * (y_c1 - y_ci);0;0];

%% Option 1 Burn directly into Transfer then into NMC

%time of flight
TOF_p = 1.834; % period # of periods determined through plug and chug
TOF = TOF_p * period; %period - sec

% Xc to obtain desired State
drift_rate = (y_ci - C*2) * TOF_p/(2*pi()); % desired drift raft to get to the target orbit within the given time - km/ period 

x_c1 = -3 * pi() /drift_rate; % change in instaneous center of ellipse

delv_1y =  [0; -n/2 * (x_c1 - x_c0);0];

delv_o1i = delv_1y; % delta V required to change the center of the ellipse, does not include plane change

% Change in zmax for both options
delv_1z = [0;0;n * (z_max1 - z_max0)];

% Final Burn to enter desired NMC for Option 1
delv_o1f = delv_1y - delv_1z + delv_s_to_NMC;

rho_dot_tgt_chase_fo1 = delv_s_to_NMC + delv_1z; %final velocity

zulu_time_o1 = zulu_time_i + TOF/3600; % zulu time at final Burn for optino 1 - time of inital Zulu + time passed - hrs

[CATS_o1,local_time1_o1] = lighting(zulu_time_o1,longitude_i,rho_tgt_chase_i,day_of_year); %CATs at second Delta V burn for option 1

%% Option 2 Burn into NMC to NMC
% travel time is half of the period
% The burn should occour at 1800 LT in order to achieve CATS = 0

wait_time_o2 = 18 - local_time_i; % Wait before inital Burn - hrs goal is to wait until 18L
travel_time_hlf_NMC = period/2; % Travel Time for half of NMC

% Delta V calcs

y_c1_o2 = (y_ci + C * 2)/2; % calculation for transition NMC centered halfway between the desired NMC and current location
delv_1x = n/2 * (y_c1_o2 - y_ci);

% There should be issue with energy matching condition because there was
% no change in relative velocity in the y direction.
%0 == 6 * n * x_c0 + 3 * delv1(2) % check to determine if there is an energy matching condition.

delv_o2i = delv_1x + delv_1z; % Inital Delta V burn for Option 2

% Change y_c1 to y_c0
delv_o2f = [-delv_1x; 0  ; 0]  + delv_s_to_NMC ; % Final Delta V burn for Option 2

% Total Delta V
DelV_total = abs(norm(delv_o2i)) + abs(norm(delv_o2f));

rho_dot_tgt_chase_fo2 = delv_s_to_NMC + delv_1z;

% CATS Check
zulu_time_o2 = zulu_time_i + wait_time_o2 + 12; % time of inital Zulu + time passed

[CATS_o2,local_time1_o2] = lighting(zulu_time_o2,longitude_i,rho_tgt_chase_f,day_of_year); %CATs at second Delta V burn for option 2

%% Output Text
time_name = "Inital Time";
time = [0 ; zulu_time_i]; % time since begining of simulation ; zulu time
name1 = "Inital Position";
units1 = "km";
name2 = "Inital Velocity";
units2 = "km/s";
name3 = "Initial CATS";
units3 = "degrees";

text_output_phase_1(name1,rho_tgt_chase_i,units1,name2,rho_dot_tgt_chase_i,units2,name3,CATS_i,units3,time_name,time)

% Delta V o1 inital
name1 = "Delta V1 for Option 1";
var1 = delv_o1i ; 
time = [0 ; zulu_time_i];
text_output_phase_1_delv(name1,var1,time)

% Delta V o1 final
name1 = "Delta V2 for Option 1";
var1 = delv_o1f ; 
time = [TOF/3600; zulu_time_o1];
text_output_phase_1_delv(name1,var1,time)

time_name = "Final Time for Option 1:";
time = [TOF/3600; zulu_time_o1]; % time since begining of simulation ; zulu time
name1 = "Final Position for Option 1";
units1 = "km";
name2 = "Final Velocity for Option 1";
units2 = "km/s";
name3 = "Final CATS for Option 1";
units3 = "degrees";

text_output_phase_1(name1,rho_tgt_chase_f,units1,name2,rho_dot_tgt_chase_fo1,units2,name3,CATS_o1,units3,time_name,time)

% Delta V o1 inital
name1 = "Delta V1 for Option 2";
var1 = delv_o2i ; 
time = [wait_time_o2 ; zulu_time_i + wait_time_o2];
text_output_phase_1_delv(name1,var1,time)

% Delta V o1 final
name1 = "Delta V2 for Option 2";
var1 = delv_o2f ; 
time = [wait_time_o2 + 12; zulu_time_o2];
text_output_phase_1_delv(name1,var1,time)

time_name = "Final Time for Option 2:";
time = [wait_time_o2 + 12; zulu_time_o2]; % time since begining of simulation ; zulu time
name1 = "Final Position for Option 2";
units1 = "km";
name2 = "Final Velocity for Option 2";
units2 = "km/s";
name3 = "Final CATS for Option 2";
units3 = "degrees";

text_output_phase_1(name1,rho_tgt_chase_f,units1,name2,rho_dot_tgt_chase_fo2,units2,name3,CATS_o2,units3,time_name,time)

%% Extrapolation

tf = [TOF; 12*3600];

tspan1= 0:60:tf(1);

Xi1d=[rho_tgt_chase_i;rho_dot_tgt_chase_i];
% create vectors for length of tspan1of HCW
Xf1 = zeros(length(tspan1),6);
rho1 = zeros(length(tspan1),1);

% Option 1
for i = 1:length(tspan1)
    Xf1(i,:) = HCW(r_tgt-6378.137,Xi1d,tspan1(i));
    rho1(i) = norm([Xf1(i,1),Xf1(i,2),Xf1(i,3)]);
end

% create vectors for length of tspan2 of HCW
tspan2 = 0:60:(tf(2));

Xi2d=[rho_tgt_chase_i,rho_dot_tgt_chase_i];

Xf2 = zeros(length(tspan2),6);
rho2 = zeros(length(tspan2),1);

%Option 2
for i = 1:length(tspan2)
    Xf2(i,:) = HCW(r_tgt-6378.137,Xi2d,tspan2(i));
    rho2(i) = norm([Xf2(i,1),Xf2(i,2),Xf2(i,3)]);
end

%% Plot

