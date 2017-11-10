% =========================================================================
% CONSTANTS
% =========================================================================
global k_a k_r k_2 k_cat K_M epsilon_0 n_0 LtVt N_CS N_AV N_0


% Adsorption Kinetics
k_a = [0.01,0.1,0.1]; % uM^-1 s^-1
k_r = [0.0058,3.8,0.05]; % s^-1

% Plasminogen Activation
k_2 = 0.3; % s^-1
K_M = 0.19; % some units

% Fibrinolysis
k_cat = 0.2; %s^-1

% Fibre parameters 
dr = 6; dth = 6; %nm
L_M =  22.5; %nm
rho_fibre = 0.28; %g/ml

% Avogadro's Number
N_AV        = 6.02E23;

% =======================================================================
% CALCULATION OF INITIAL CLOT PROPERTIES
% =======================================================================

% Detrmining number of cross-sections
cord = dr;
constraint = 0;
it_check=0;
group = 1;
while constraint < rho_0/1000
    % Update the iteration 
    it_check = it_check + 1;
    % Assign random value 
    %R_f0(it_check) = round(floor(normrnd(mu,sigma))/cord)*cord;
    R_f0(it_check) = round(normrnd(mu,sigma));
    
    % Update constraint
    consStore = constraint;
        
    constraint = group*sum(R_f0.^2)*L_M*rho_fibre*pi()/(U_0*1E9);
   
    
    % Break the loop in case of error
    if it_check > 100000
       
        error('ERROR: Calculated number of cross-sections is too big. Not enough memory');
        
    end
end

N_CS = it_check;

% Keeping the constraint of the fibrin density
Rf_end = sqrt((rho_0/1000-consStore)*U_0*1E9/(rho_fibre*pi()*L_M*group));
R_f0(end) = round(Rf_end/cord)*cord;

% Determining binding locations per cross-section
N_0 = zeros(1,N_CS);

for it_out = 1:N_CS

    for it = 1: round(R_f0(it_out)/dr)
        
        N_0(it_out) = N_0(it_out) + pi()/asin(dth/(2*it*dr));
        
    end

end

% Determining initial binding site concentration
n_0 = group*sum(N_0)*1E6*1E15/(U_0*N_AV);


% % Fibre Length Density
% LtVt = (rho_0/(rho_fibre*1000))*10^21/(pi()*R_f0.^2); %nm/cm^3
% 
% % Cross-sections per volume 
% CS_V = LtVt/(1E12)/L_M; % um^-3
% 
% 
% % Finding the fibrin volume fraction
% phi_f = pi()*R_f0^2*LtVt*1E-21;

phi_f = rho_0/rho_fibre;

% Finding the porosity
epsilon_0 = 1 - phi_f; %unitless


