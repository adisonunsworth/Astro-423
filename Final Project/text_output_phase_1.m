function text_output_phase_1(name1,var1,units1,name2,var2,units2,name3,var3,units3,time_name,time)

% name1   = "Inital Position";
% var1    = "var1";
% time   = [1,15];
% name2   = "name2";
% var2    = "var2";
% 
% name3   = "name3";
% var3    = "var3";

%TEXT_OUTPUT Summary of this function goes here
%   Detailed explanation goes here
fprintf('At %s %d [hrs] and zulu time %d [hrs]: \n',time_name,time(1),time(2))
fprintf('   %s in the RIC frame is:\n   %5d x-hat [%s]\n   %5d y-hat [%s]\n   %5d z-hat [%s]\n\n',name1,var1(1),units1,var1(2),units1,var1(3),units1)

fprintf('   %s in the RIC frame is:\n   %5d x-hat [%s]\n   %5d y-hat [%s]\n   %5d z-hat [%s]\n\n',name2,var2(1),units2,var2(2),units2,var2(3),units2)

fprintf('   %s  is:\n   %5.2d x-hat [%s]\n\n',name3,var3,units3)



% fprintf('Position in the RIC frame after [s]:\n %d sec\n [km]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(1), Xf2(1), Xf2(2), Xf2(3))
% fprintf('Velocity in the RIC frame after [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(1),Xf2(4), Xf2(5), Xf2(6))
% 
% fprintf('Delta V_2 in the RIC frame at [s]:\n %d sec\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(1) , dV2(1), dV2(2), dV2(3))
% fprintf('Position in the RIC frame after [s]:\n %d sec\n [km]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(2), Xf3(1), Xf3(2), Xf3(3))
% fprintf('Velocity in the RIC frame after [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(2),Xf3(4), Xf3(5), Xf3(6))
% fprintf('Time before final burn [s]:\n   %d sec\n\n',tf(2) - tf(1))
% 
% fprintf('Delta V_3 in the RIC frame at at [s]:\n %d sec\n [km/s]:\n   %d x-hat\n   %d y-hat\n   %d z-hat\n\n',tf(2), -Xf3(4), -Xf3(5), -Xf3(6))


%end

