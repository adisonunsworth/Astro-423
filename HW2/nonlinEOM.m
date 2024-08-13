function Xdot = nonlinEOM(~,X)

mu=398600.5; %Gravitational parameter for Earth

% The state vector x is defined as follows:
x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);
r_tgt = X(7);
r_tgtdot = X(8);
theta_tgt = X(9);
theta_tgtdot= X(10);

Xdot = [0 0 0 0 0 0 0 0 0 0]'; % Initialize State Vector derivative

%Define state vector derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - You complete the code below. 
%   - Add additional lines before xdot(1) to define terms that are common 
%   to multiple state vector component derivaties if desired.
%   - All state vector component derivative equations must be written in
%   terms of ONLY state vector components (i.e., x(1), x(2), x(3), etc) and
%   mu.


A = ((r_tgt+x)^2 + y^2 + z^2)^(3/2);
theta_tgt_doubledot = -2 * r_tgtdot * theta_tgtdot/ r_tgt;

Xdot(1)= xdot;   % xdot
Xdot(2)= ydot;   % ydot
Xdot(3)= zdot;   % zdot
Xdot(4)= -mu *(r_tgt + x) / (A) + mu / (r_tgt^2) + 2 * theta_tgtdot * ydot + theta_tgt_doubledot * y + theta_tgtdot^2 * x;  % x doubledot
Xdot(5)= -mu * y / A - 2 * theta_tgtdot * xdot - theta_tgt_doubledot * x + theta_tgtdot^2 * y; % y doubledot
Xdot(6)= -(mu * z)/A ;  % z doubledot
Xdot(7)= X(8);   %r_tgt dot
Xdot(8)= r_tgt * theta_tgtdot^2 - mu / r_tgt^2;   % r_tgt double dot
Xdot(9)= theta_tgtdot;   % theta_tgt dot
Xdot(10)= theta_tgt_doubledot;  % theta_tgt doubledot

Xdot = double(Xdot); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end