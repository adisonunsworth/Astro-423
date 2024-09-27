function text_output_phase_1(name1,var1,units1,name2,var2,units2,name3,var3,units3,time_name,time)
%TEXT_OUTPUT Summary of this function goes here
%   Detailed explanation goes here
fprintf('At %s %d [hrs] and zulu time %d [hrs]: \n',time_name,time(1),time(2))
fprintf('   %s in the RIC frame is:\n   %5d x-hat [%s]\n   %5d y-hat [%s]\n   %5d z-hat [%s]\n\n',name1,var1(1),units1,var1(2),units1,var1(3),units1)

fprintf('   %s in the RIC frame is:\n   %5d x-hat [%s]\n   %5d y-hat [%s]\n   %5d z-hat [%s]\n\n',name2,var2(1),units2,var2(2),units2,var2(3),units2)

fprintf('   %s  is:\n   %5.2d x-hat [%s]\n\n',name3,var3,units3)
end

