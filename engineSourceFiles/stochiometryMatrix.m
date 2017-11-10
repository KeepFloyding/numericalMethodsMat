% ======================================================================
% STOCHIOMETRY MATRIX
% ======================================================================

v = zeros(N_CS*N_rxns + 1,N_CS*N_species + 1);

for s_it = 1:N_CS
    
    v(s_it + 0*N_CS,s_it + 2*N_CS) = -1; % PLS desorption
    v(s_it + 1*N_CS,s_it + 2*N_CS) =  1; % PLS adsorption
    
    v(s_it + 2*N_CS,s_it + 2*N_CS)   = 1;  % M-M rxn
    v(s_it + 2*N_CS,s_it + 1*N_CS)   = -1;
    
    v(s_it + 3*N_CS,s_it + 3*N_CS)   = -1; % Cat rxn
    v(s_it + 3*N_CS,s_it + 2*N_CS)   = -1;
    v(s_it + 3*N_CS,s_it + 4*N_CS)   =  1;
    
    v(s_it + 4*N_CS,s_it + 4*N_CS)   = -1; % L_PLS dis
    
    v(s_it + 5*N_CS,s_it) = 1;             % tPA adsorption
    v(s_it + 6*N_CS,s_it) = -1;            % tPA desorption
    
    v(s_it + 7*N_CS,s_it + 1*N_CS) =  1; % PLG adsorption
    v(s_it + 8*N_CS,s_it + 1*N_CS) = -1; % PLG desorption
    
   
end

%                       Stochiometry Matrix
%  rxn/species    n_tPA   n_PLG   n_PLS   n_tot   L_PLS
%  0   PLS des     0       0       -1      0       0
%  1   PLS ads     0       0       +1      0       0
%  2   M-M rxn     0       -1      +1      0       0
%  3   Cat rxn     0       0       -1      -1      +1
%  4   L_PLS des   0       0       0       0       -1
%  5   tPA ads     +1      0       0       0       0
%  6   tPA des     -1      0       0       0       0
%  7   PLG ads     0       +1      0       0       0
%  8   PLG des     0       -1      0       0       0
    
