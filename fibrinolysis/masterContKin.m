% Continuous simulation of the kinetic equations, lysis is defined as a
% homogeneous degradation of the fibre rather than discrete cutting
clear
clc


% =======================================================================
% INPUTS
% =======================================================================
global R_f0

% Kinetics
R_f0 = 200;
rho_0 = 0.5; % mg/ml
C_tPA0 = 0.04; C_PLG0 = 2.2; C_PLS0 = 0;

% Simulation
startTime = 0; %s 
finalTime = 10*60; % s
nt = 1000;

% =======================================================================
% CONSTANTS
% =======================================================================

global k_a k_r k_2 K_M k_cat epsilon_0 n_0 LtVt

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

% Avogadro's constant
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
epsilon_0 = 1 - 3/280; %unitless

% =======================================================================
% ODE SOLVER
% =======================================================================

timespan = linspace(startTime, finalTime, nt); % Discretised timeline

% Specifying initial condtions
% 1: C_tPA; 2: n_tPA; 3: C_PLG; 4: n_PLG; 5: C_PLS; 6: n_PLS; 7: n_tot; 8:
% L_PLS

IC = zeros(8,1);
IC(1) = C_tPA0*epsilon_0;
IC(3) = C_PLG0*epsilon_0;
IC(5) = C_PLS0*epsilon_0;
IC(7) = n_0*(1-epsilon_0);

[t, sol] = ode45(@odeEquations, timespan, IC);

% Outputing variables
C_tPA = sol(:,1);
n_tPA = sol(:,2);
C_PLG = sol(:,3);
n_PLG = sol(:,4);
C_PLS = sol(:,5);
n_PLS = sol(:,6);
n_tot = sol(:,7);
L_PLS = sol(:,8);

% Calculating corresponding properties
epsilon =  clotProperties(n_tot);

% plot(t,epsilon)
% hold on
% =======================================================================
% POST-PROCESS
% =======================================================================
% M-file that does all the plotting. 
plotFile;
hold all



