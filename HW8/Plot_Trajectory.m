function  [] = Plot_Trajectory(r_tgt,Xi_HCW,tspan)

mu = 398600.5;


v_tgt   = sqrt(mu/r_tgt);

    v_chase = sqrt(2 *(mu/(r_tgt+x) - mu/(2*r_tgt)));
    
    ydot_I =  v_chase - v_tgt;

    n = sqrt(mu/r_tgt^3);
    
    y = 0 ; %initial in-track position
    z = 0 ; %initial cross-track position
    
    xdot_ijk = 0; %initial radial velocity (wrt the IJK resolved in the RIC frame)
    ydot_ijk = ydot_I; %initial in-track velocity (wrt the IJK resolved in the RIC frame)
    zdot_ijk = 0; %initial cross-track velocity (wrt the IJK resolved in the RIC frame)
    


    V = [xdot_ijk; ydot_ijk; zdot_ijk] + cross([0;0;-n], [x;0;0]);
    
    xdot = V(1);
    ydot = V(2);
    zdot = V(3);




counter = 1;

Xi_HCW = 

X_HCW=inf(length(tspan),6); %Pre-allocate matrix to store HCW trajectory
    for ii=1:length(tspan)
        X_HCW(counter,:) = HCW(r_tgt-6378.137,Xi_HCW,tspan(ii));  %%%% YOU WRITE THE FUNCTION HCW.m
        counter = counter+1;
    end

figure
    plot(X_HCW(:,2),X_HCW(:,1),'--')
    xlabel('y - in-track [km]')
    ylabel('x - radial [km]')
    legend('Numerically Integrated Trajectory (ode45)','HCW Trajectory','Location','southwest')
    grid on