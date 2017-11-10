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

% Fibre Length Density
LtVt = (rho_0/(rho_fibre*1000))*10^21/(pi()*R_f0^2); %nm/cm^3

% Cross-sections per volume 
CS_V = LtVt/(1E12)/L_M; % um^-3

% Binding site per cross-section of size R_f0
BS_CS = 0;
for it = 1: round(R_f0/dr)
    
    BS_CS = BS_CS + pi()/asin(dth/(2*it*dr));
    
end

% Total concentration of binding sites
n_0 = CS_V * BS_CS*1E6*1E15/N_AV; % uM

% Finding the fibrin volume fraction
phi_f = pi()*R_f0^2*LtVt*1E-21;

% Finding the porosity
epsilon_0 = 1 - phi_f; %unitless

% Total number of cross-sections
N_CS = round(CS_V*U_0);

% Total number of binding sites
N_0 = CS_V*BS_CS*U_0;
