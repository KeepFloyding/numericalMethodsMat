% Runge-Kutta 4th order solver, no edits need to be made to this file. 
% Differential equations should be located in file known as gradient.

function sol = rungeKuttaFourth(y, x, dx)

solve = @gradient;

k_1 = dx*solve(x,y);
k_2 = dx*solve(x + dx/2,y + k_1/2);
k_3 = dx*solve(x + dx/2,y + k_2/2);
k_4 = dx*solve(x + dx, y + k_3);


sol = y + k_1/3 + k_2/6 + k_3/6 + k_4/3;


end