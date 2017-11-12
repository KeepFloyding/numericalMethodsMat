function eqns = odeEquations(t,y)

eqns = zeros(8,1);

global k_a k_r k_2 K_M k_cat 

% =======================================================================
% EQUATION ORDER
% =======================================================================

Cs_tPA = y(1); % 1: Cs_tPA
ns_tPA = y(2); % 2: ns_tPA
Cs_PLG = y(3); % 3: Cs_PLG
ns_PLG = y(4); % 4: ns_PLG
Cs_PLS = y(5); % 5: Cs_PLS
ns_PLS = y(6); % 6: ns_PLS
ns_tot = y(7); % 7: ns_tot
L_PLS =  y(8); % 8: L_PLS

% =======================================================================
% AUXILLARY FUNCTIONS
% =======================================================================

epsilon =  clotProperties(ns_tot);

ns_free = ns_tot - (ns_tPA + ns_PLG + ns_PLS);

% =======================================================================
% EQUATIONS
% =======================================================================

eqns(1) = -(k_a(1)*(Cs_tPA/epsilon)*ns_free - k_r(1)*ns_tPA) ;
eqns(2) =  (k_a(1)*(Cs_tPA/epsilon)*ns_free - k_r(1)*ns_tPA);
eqns(3) = -(k_a(2)*(Cs_PLG/epsilon)*ns_free - k_r(2)*ns_PLG);
eqns(4) =   k_a(2)*(Cs_PLG/epsilon)*ns_free - k_r(2)*ns_PLG - k_2*ns_PLG*ns_tPA/(K_M*(1-epsilon) + ns_PLG);
eqns(5) = -(k_a(3)*(Cs_PLS/epsilon)*ns_free - k_r(3)*ns_PLS) + k_r(3)*L_PLS*(1-epsilon);
eqns(6) =   k_a(3)*(Cs_PLS/epsilon)*ns_free - k_r(3)*ns_PLS - k_cat*ns_PLS + k_2*ns_PLG*ns_tPA/(K_M*(1-epsilon) + ns_PLG);
eqns(7) =  -k_cat*ns_PLS;
eqns(8) = - k_r(3)*L_PLS*(1-epsilon) + k_cat*ns_PLS;

% =======================================================================
% CONDITIONALS TO AVOID ASYMPTOTES
% =======================================================================

if (epsilon == 1)
   
    eqns(4) =   k_a(2)*(Cs_PLG/epsilon)*ns_free - k_r(2)*ns_PLG;
    eqns(6) =   k_a(3)*(Cs_PLS/epsilon)*ns_free - k_r(3)*ns_PLS - k_cat*ns_PLS;
    
    
end


end