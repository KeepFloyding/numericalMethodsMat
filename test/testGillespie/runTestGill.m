% Test script for Gillespie solver
%   Equation is the reversible reaction A + B <-> C
%   Forward rate constant k_f = 1; mol^-1 s^-1
%   Reverse rate constant k_r = 1; s^-1
%   Intial concentration of A, B and B is 1, 0.5 and 0 moles
%   propensityFunction.m contains the kinetic equations

% Simulation settings
delta_t = 0.01; % Units are seconds
nt = 100; 

% Storage
t = zeros(nt);
C_a = zeros(nt);
C_b = zeros(nt);
C_c = zeros(nt);

% Global variables
global k_f k_r

% Factor to turn moles into molecules
U_0 = 1; %volume of system;
N_AV = 6.02E23; % Avogadro's number
F = U_0*1E-6*1E-15*N_AV;

k_f = 1/F;
k_r = 1;

% Stochiometry Matrix
v = zeros(2,3);

%  rxn/species      A         B         C 
%  0   Forward      -1        -1       +1  
%  1   Backward     +1        +1       -1  

v(1,:) = [-1,-1,1];
v(2,:) = [1,1,-1];

% State arrray (initial conditions)
X = [1*F, 0.5*F, 0];
C_a(1) = 1;
C_b(1) = 0.5;
C_c(1) = 0;

for it = 2:nt
    
    % Update the current time
    current_t = (it-1)*delta_t;
    t(it) = current_t;
    
    % Run the Gillespie solver
	X_nx = gillespieSolver(X,v, delta_t,@propensityFunction);
    
    % Store variables
    C_a(it) = X_nx(1)/F;
    C_b(it) = X_nx(2)/F;
    C_c(it) = X_nx(3)/F;
    
    % Update the state array
	X = X_nx;
    
end

% Plot the results
plot(t,C_a,t,C_b,t,C_c)
legend('A','B','C')
