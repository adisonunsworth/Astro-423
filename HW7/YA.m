function [Xtraj]=YA(X0,nu0,e,a,tspan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  YA.m
%%
%%  This function propagates an initial scaled relative motion state for 
%%  the chase forward in time over a specified timespan using the Yamanaka 
%%  and Ankersen state transition matrix to generate the chase's relative 
%%  motion trajectory in scaled coordinates. 
%%
%%  Author        : Lt Col Greg Frey                   Sep 2023
%%  Modified by   :
%%
%%  [Xtraj]=YA(X0,nu0,e,a,tspan)
%%
%%  Inputs        :
%%    X0          : Scaled Initial state for the chase,
%%                  [xbar; xbar'; ybar; ybar'; zbar; zbar'], 
%%                  units of km and km/rad
%%    nu0         : True anomaly of the target at t=0 [rad]
%%    e           : Target orbit eccentricity [-]
%%    a           : Target orbit semi-major axis [km]
%%    tspan       : Vector of times used to calculate trajectory,
%%                  [0:timestep:final time] [sec]
%%
%%  Outputs       :
%%    X_traj      : Matrix containing Relative motion trajectory and
%%                  associated data. 
%%                  Each row contains [xbar xbar' ybar ybar' zbar zbar' t nu]
%%
%%  Coupling
%%   newton.m    Calculates eccentric anomaly of the target at each timestep
%%
%%  References             :
%%    Yamanaka, Koji, and Finn Ankersen. "New state transition matrix for 
%%    relative motion on an arbitrary elliptical orbit." J Guid Control Dyn 
%%    25.1 (2002): 60-66.
%%    Spacecraft Formation Flying by Alfriend et. al., Section 5.6
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize array to store the trajectory and data
Xtraj=inf(length(tspan),8);

%Define constants
mu=398600.5;    %Gravitational parameter for Earth

% %Calculate additional orbital parameters for the target
n=sqrt(mu/a^3);             %Target mean motion [rad/s]
p=a*(1-e^2);                %Parameter for the target [km]
h=sqrt(mu*p);               %Specific angular momentum of target [km^2/s]

%Calculate the relative motion trajectory using the Yamanaka and Ankersen STM
for ii=1:length(tspan)
    t=tspan(ii); %Extract the current time
    
    % Calculate initial value for Mean anomaly at t=0 
    E0=acos((e+cos(nu0))/(1+e*cos(nu0)));   %Eccentric anomaly [rad]
    % Half-plane check for E
        if nu0>pi
             E0 = 2*pi-E0;
        end
    M0=E0-e*sin(E0);   
    
    % Determine target true anomaly at the current time
    M=M0+n*t;          %Calculate mean anomaly [rad]
    M=mod(M,2*pi);     %Revcheck M to ensure between 0 and 2pi
    [E] = newton(M,e); %Solve Kepler's problem to obtain Eccentric anomaly using Newton's method [rad]
    nu = acos((cos(E)-e)/(1-e*cos(E))); %Calculate true anomaly [rad]
    % Half-plane check
    if E>pi
         nu = 2*pi-nu;
    end
    
    %Propagate state forward in time using the YA STM
    %Calculate terms to simplify STM
    I=(mu^2/h^3)*t;
    eta=sqrt(1-e^2);
    
    kt=1+e*cos(nu);                         %parameter k evalulated at the current nu value
    st=kt*sin(nu);                          %parameter s evalulated at the current nu value
    stp=cos(nu)+e*(cos(nu)^2-sin(nu)^2);    %derivative of s wrt nu, evalulated at the current nu value
    ct=kt*cos(nu);                          %parameter c evalulated at the current nu value
    ctp=-sin(nu)-2*e*sin(nu)*cos(nu);       %derivative of c wrt nu, evalulated at the current nu value
    
    k0=1+e*cos(nu0);                        %parameter k evalulated at the initial nu value (nu0)
    s0=k0*sin(nu0);                         %parameter s evalulated at the initial nu value (nu0)
    c0=k0*cos(nu0);                         %parameter c evalulated at the initial nu value (nu0)
    
    %Form matricies used to calculate the STM
    phi_nu=[st             ct             2-3*e*st*I            0    0        0;
            stp            ctp           -3*e*(stp*I+(st/kt^2)) 0    0        0;
            ct*(1+(1/kt)) -st*(1+(1/kt)) -3*kt^2*I              1    0        0;
            -2*st          e-2*ct        -3*(1-2*e*st*I)        0    0        0;
            0              0              0                     0    cos(nu)  sin(nu);
            0              0              0                     0   -sin(nu)  cos(nu)];
        
    iphi_0=(1/eta^2)*[-3*s0*((k0+e^2)/k0^2)    c0-2*e  0      -s0*((k0+1)/k0)       0                0;
                      -3*(e+c0/k0)            -s0      0      -(c0*((k0+1)/k0)+e)   0                0;
                       3*k0-eta^2              e*s0    0       k0^2                 0                0;
                      -3*e*s0*((k0+1)/k0^2)   -2+e*c0  eta^2  -e*s0*((k0+1)/k0)     0                0;
                       0                       0       0       0                    eta^2*cos(nu0)  -eta^2*sin(nu0);
                       0                       0       0       0                    eta^2*sin(nu0)   eta^2*cos(nu0)];
    
    STM=phi_nu*iphi_0;   %YA State Transition Matrix (STM)
    Xt=STM*X0; %Calulate the chase's state at the current time. 
    
    %Store the chase's state vector, the current time and the target's true 
    %anomaly in the trajectory array
    Xtraj(ii,:)=[Xt' t nu];
 end