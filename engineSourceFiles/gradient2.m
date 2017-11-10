function eqns = gradient2(x,y)

% ========================================================================
% CALLING CONSTANTS AND OTHER REQUIRED INFORMATION
% ========================================================================

global k_a k_r store extentLysis alpha

ns_tPA = store.ns_tPA;
ns_PLG = store.ns_PLG;
ns_PLS = store.ns_PLS;
ns_tot = store.ns_tot;
L_PLS = store.L_PLS;
E_L = extentLysis;

% ========================================================================
% INITIALISING EQUATIONS AND VARIABLES
% ========================================================================

eqns = zeros(4,1);
Cs_tPA = y(1);
Cs_PLG = y(2);
Cs_PLS = y(3);
cell = y(4);

if E_L== 1
    E_L = 0.95;
end

K_rel = alpha*E_L/(1 - E_L);

if cell < 0
    cell = 0;
    
end

% ========================================================================
% CONDITIONALS
% ========================================================================

if (ns_PLG < 0) 
    ns_PLG = 0;
   
elseif (ns_tPA < 0)
    ns_tPA = 0;
    
elseif (ns_PLS < 0)
    ns_PLS = 0;
    
end


 ns_free = ns_tot - ns_PLG - ns_tPA - ns_PLS;

 if ns_free < 0
     
     ns_free = 0;
     
 end
 
% ========================================================================
% DIFFERENTIAL EQUATIONS
% ========================================================================
    
eqns(1) = -(k_a(1)*(Cs_tPA)*ns_free - k_r(1)*ns_tPA);
eqns(2) = -(k_a(2)*(Cs_PLG)*ns_free - k_r(2)*ns_PLG);
eqns(3) = -(k_a(3)*(Cs_PLS)*ns_free - k_r(3)*ns_PLS) ...
     + k_r(3)*L_PLS;
eqns(4) = -K_rel*cell;

end