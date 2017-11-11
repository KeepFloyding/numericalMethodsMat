function eqns = gradient(x,y)

% -------------------------------------------------------------------------
% DIFFERENTIAL EQUATIONS
% -------------------------------------------------------------------------
 
global k_f k_r

eqns = zeros(3,1);

C_a = y(1);
C_b = y(2);
C_c = y(3);

%   Equation is the reversible reaction A + B <-> C

eqns(1) = -k_f*C_a*C_b + k_r*C_c;
eqns(2) = -k_f*C_a*C_b + k_r*C_c;
eqns(3) =  k_f*C_a*C_b - k_r*C_c;

end