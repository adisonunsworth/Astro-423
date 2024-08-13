%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A423_HW8.m
%
%  This script does what is required for HW8 
%
%  Author        :  C2C Unsworth                           2024
% 
%  Inputs        :
%
%   Outputs       :
%       X_traj     : Relative motion trajectory [km, km/rad, s, and rad].
%
%  Coupling
%           HCW_Targeting_without_dV2
%   
%
%  References             :
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
format long

%% Target vehicle information
r_tgt=7000; %Geo orbit radius [km]

%% ICs for the Chase
% Initial position [km]
x0= 0;
y0= -5;
z0= 0;
r0=[x0; y0; z0];
% initial velocity km/s
xd0= 0;
yd0= 0;
zd0= 0;
v0=[xd0; yd0; zd0];
% Pre DelV Initial state vector
X01=[r0; v0];

%% Fill in the following information about the desired trajectory 
TOF1 = 3600*5; % Time from now until reaching desired location [seconds]
TOF2 = 3600*12;

tspan1= 0:60:TOF1; % Total Time operation should take
tspansim1 = 0:60:TOF2; % Total Simulation Time

%% Setup for 
% Desired final position (km)
xf = [1,0];
yf = [-3.5,-1];
zf = [0,0];
rf = [xf;yf;zf];

% Desired final velocity (km/s)
xdf= 0;
ydf= 0;
zdf= 0;
vf2=[xdf; ydf; zdf];

% Desired final time for each burn
tf = [1*3600,5*3600]; %inital burn time and final burn time

%Desired final position vector
Xf1d = [rf(:,1)];

Xf2d = [rf(:,2)];

% Find Inital DelV and final State Vector
[dV1, Xf02]= HCW_Targeting_without_dV2(r_tgt, X01, Xf1d, tf(1));

% desired initial velocity
xd0= dV1(1);
yd0= dV1(2);
zd0= dV1(3);
v0=[xd0; yd0; zd0];

tspan1= 0:60:tf(1);

Xi1d=[r0,v0];
% create vectors for length of tspan1of HCW
Xf1 = zeros(length(tspan1),6);
rho1 = zeros(length(tspan1),1);

% finds state vectors for chase depending on inital Conditions
for i = 1:length(tspan1)
    Xf1(i,:) = HCW(r_tgt-6378.137,Xi1d,tspan1(i));
    rho1(i) = norm([Xf1(i,1),Xf1(i,2),Xf1(i,3)]);
end
% Creates vectors for simulated trajectory
Xf1sim = [Xf1;zeros(length(tspansim1) - length(tspan1),6)];
rho1sim = [rho1;zeros(length(tspansim1) - length(rho1),1)];

for k = tf(1)/60:length(tspansim1) % uses the first time to determine when to start the simulation
    Xf1sim(k,:) = HCW(r_tgt-6378.137,Xi1d,tspansim1(k));
    rho1sim(k) = norm([Xf1sim(k,1),Xf1sim(k,2),Xf1sim(k,3)]);
end

%Plots
figure(1)
title("Inital Delta V Burn with Simulation to 12 hours")
hold on
p1_3_sim = plot3(Xf1sim(:,2),Xf1sim(:,1),Xf1sim(:,3),'LineStyle','--','Color','blue'); % plots simulated plot
p1_3 = plot3(Xf1(:,2),Xf1(:,1),Xf1(:,3));
hold off
grid on
ylabel("distance x axis (km)")
xlabel("distance y axis (km)")
zlabel("distance z axis (km)")
legend({'Simulation to 12 hrs','Simulation to Second Delta V'},'Location','southeast')
legend('boxoff')

figure(2)
hold on
title("Inital Delta V Burn with Simulation to 12 hours")
p1_rho_sim = plot(tspansim1,rho1sim,'LineStyle','--','Color','blue');
p1_rho = plot(tspan1,rho1);
hold off
grid on
ylabel("distance between chase and tgt km")
xlabel("time s")
legend({'Simulation to 12 hrs','Simulation to Second Delta V'},'Location','southeast')
legend('boxoff')

minimum_dis = min(rho1sim);
txt = ['Minimum Distance: ' num2str(minimum_dis) ' km'];
text(500,10,txt)



%%
[dV2, Xf3]=HCW_Targeting_without_dV2(r_tgt, Xf02, Xf2d, tf(2) - tf(1)); 

% desired initial velocity

% Xf2 Inital Position of the Second burn
x = Xf02(1);
y = Xf02(2);
z = Xf02(3);
r2 = [x;y;z];

%velocity after second burn
% xd0= dV2(1);
% yd0= dV2(2);
% zd0= dV2(3);
xd0= Xf02(4) + dV2(1);
yd0= Xf02(5) + dV2(2);
zd0= Xf02(6) + dV2(3);
v2 = [xd0; yd0; zd0];

tspan2 = 0:60:(tf(2)-tf(1));
tspansim2 = 0:60:(TOF2 -tf(1) );

Xi2d=[r2,v2];
% create vectors for length of tspan2 of HCW
Xf2 = zeros(length(tspan2),6);
rho2 = zeros(length(tspan2),1);

for i = 1:length(tspan2)
    Xf2(i,:) = HCW(r_tgt-6378.137,Xi2d,tspan2(i));
    rho2(i) = norm([Xf2(i,1),Xf2(i,2),Xf2(i,3)]);
end

% Creates vectors for simulated trajectory
Xf2sim = [Xf2;zeros(length(tspansim2) - length(tspan2),6)]; 
rho2sim = [rho2;zeros(length(tspansim2) - length(rho2),1)];
% Runs through Simulated Trajectory
for k = (tf(2) - tf(1))/60:length(tspansim2) % uses the first time to determine when to start the simulation
    Xf2sim(k,:) = HCW(r_tgt-6378.137,Xi2d,tspansim2(k));
    rho2sim(k) = norm([Xf2sim(k,1),Xf2sim(k,2),Xf2sim(k,3)]);
end

figure(3)
title("Second Delta V Burn with with Simulation to 12 hours")
hold on
p2_3_sim = plot3(Xf2sim(:,2),Xf2sim(:,1),Xf2sim(:,3),'LineStyle','--','Color','blue'); % plots simulated plot
p2_3 = plot3(Xf2(:,2),Xf2(:,1),Xf2(:,3));
hold off
grid on
ylabel("distance x axis (km)")
xlabel("distance y axis (km)")
zlabel("distance z axis (km)")
legend({'Simulation to 12 hrs','Simulation to Third Delta V'},'Location','southeast')
legend('boxoff')

figure(4)
hold on
title("Second Delta V Burn with Simulation to 12 hours")
p2_rho_sim = plot(tspansim2,rho2sim,'LineStyle','--','Color','blue');
p2_rho = plot(tspan2,rho2);
hold off
grid on
ylabel("Distance Between Chase and Tgt (km)")
xlabel("Time (s)")
legend({'Simulation to 12 hrs','Simulation to Third Delta V'},'Location','southeast')
legend('boxoff')

minimum_dis = min(rho2sim);
txt = ['Minimum Distance: ' num2str(minimum_dis) ' km'];
text(500,30,txt)

%% Prints Relevant Answers

fprintf('Delta V_1 in the RIC frame at [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',0, dV1(1), dV1(2), dV1(3))
fprintf('Position in the RIC frame after [s]:\n %d sec\n [km]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(1), Xf2(1), Xf2(2), Xf2(3))
fprintf('Velocity in the RIC frame after [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(1),Xf2(4), Xf2(5), Xf2(6))

fprintf('Delta V_2 in the RIC frame at [s]:\n %d sec\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(1) , dV2(1), dV2(2), dV2(3))
fprintf('Position in the RIC frame after [s]:\n %d sec\n [km]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(2), Xf3(1), Xf3(2), Xf3(3))
fprintf('Velocity in the RIC frame after [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(2),Xf3(4), Xf3(5), Xf3(6))
fprintf('Time before final burn [s]:\n   %d sec\n\n',tf(2) - tf(1))

fprintf('Delta V_3 in the RIC frame at at [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(2), -Xf3(4), -Xf3(5), -Xf3(6))

%% Combined Plot
% 
figure(5)
hold on
% plot3([Xf1(:,2);Xf2(:,2)],[Xf1(:,1);Xf2(:,1)],[Xf1(:,3);Xf2(:,3)],'blue') % p1_3;

plot3(Xf1(:,2),Xf1(:,1),Xf1(:,3),'blue')
plot3(Xf2(:,2),Xf2(:,1),Xf2(:,3),'black')

title('Combined Trajectory')


ylabel("distance x axis (km)")
xlabel("distance y axis (km)")
zlabel("distance z axis (km)")
Target = plot3(0,0,0,'r*');

legend({'First Trajectory','Second Trajectory'},'Location','southeast')
legend('boxoff')

hold off

% 
% figure(6)
% hold on
% plot(tspan1,rho1) % p1_rho;
% plot(tspan2,rho2) % p2_rho;
% hold off

% I have spent far too much time on this HW, almost as much as a
% propulsion project, I hope my graphs are good enough without a combined
% plot.