% Test script for Runge-Kutta solver
%   Equation is the reversible reaction A + B <-> C
%   Forward rate constant k_f = 1; mol^-1 s^-1
%   Reverse rate constant k_r = 1; s^-1
%   Intial concentration of A, B and B is 1, 0.5 and 0 moles
%   gradient.m contains the kinetic equations

% INPUTS

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
k_f = 1;
k_r = 1;

% Initial Conditions
C_a(1) = 1;
C_b(1) = 0.5;
C_c(1) = 0;

input = [C_a(1);C_b(1);C_c(1)];

for it = 2:nt
    
    % Update the current time
    current_t = (it-1)*delta_t;
    t(it) = current_t;
    
    % Run Runge Kutta Solver
    sol = rungeKuttaFourth(input, current_t, delta_t);

    % Store variables
    C_a(it) = sol(1);
    C_b(it) = sol(2);
    C_c(it) = sol(3);

    % Reset solver input
	input = [sol(1); sol(2); sol(3)];

end
 
% Plot the results
plot(t,C_a,t,C_b,t,C_c)
legend('A','B','C')

