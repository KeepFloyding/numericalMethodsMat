function epsilon =  clotProperties(ns_tot)

    global R_f0 n_0 epsilon_0 LtVt
            
    % Radius of the fibrin fibre
    
    %R_f = real(R_f0*sqrt(ns_tot/(n_0*(1-epsilon_0))));
    
    
    % Finding the porosity
    %epsilon = 1 - pi()*R_f.^2*(LtVt)*10^-21; %unitless

    EL = (1 - ns_tot/(n_0*(1-epsilon_0)));
    
    phi_f = (1-epsilon_0)*(1-EL);
    
    epsilon = 1 - phi_f;
    
    
end
