function B=axisrot(A,axis,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Use           : B=axisrot(A,axis,alpha)
%%
%% This function performs a rotation of angle ALPHA about a desired axis.
%%
%%   Author       : Dr. RON LISOWSKI, DFAS,      5 Jan 95
%%   In MatLab    : Thomas L. Yoder, LtC, USAFA, Spring 00
%%
%%   Input        :
%%     A          % Input vector                            Vector of dimension three
%%     axis       % desired axis for rotation:              1, 2 or 3
%%     alpha      % Angle of rotation                       radians
%%
%%   Output       :
%%     B          % Rotated Vector                          Vector of dimension three
%%
%%   Locals       : None.
%%
%%   Coupling     :
%%     mag        % Finds the magnitude of a vector
%%
%%   References   :
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch axis
case 1   
   %rotate about the 1st axis
   B(1)=A(1);
   B(2)=A(2)*cos(alpha)+A(3)*sin(alpha);
   B(3)=-A(2)*sin(alpha)+A(3)*cos(alpha);
case 2
   % rotate about the 2nd axis
   B(1)=A(1)*cos(alpha)-A(3)*sin(alpha);
   B(2)=A(2);
   B(3)=A(1)*sin(alpha)+A(3)*cos(alpha);
case 3
   % rotate about the 3rd axis
   B(1)=A(1)*cos(alpha)+A(2)*sin(alpha);
   B(2)=-A(1)*sin(alpha)+A(2)*cos(alpha);
   B(3)=A(3);
otherwise
   disp('AxisRot axis number not 1, 2 or 3')
   B = A;
end
   % Ensure B is a column vector
   B = B(:);