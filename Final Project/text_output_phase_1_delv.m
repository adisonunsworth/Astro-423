function text_output_phase_1_delv(name1,var1,time)
%TEXT_OUTPUT Summary of this function goes here
%   Detailed explanation goes here
fprintf('%s occours at %d [hrs] and zulu time %d [hrs]: \n',name1,time(1),time(2))
fprintf('   %s in the RIC frame is:\n   %5d x-hat [km/s]\n   %5d y-hat [km/s]\n   %5d z-hat [km/s]\n\n',name1,var1(1),var1(2),var1(3))
end

