% ======================================================================
% STOCHIOMETRY MATRIX
% ======================================================================

v = zeros(N_rxns,N_species);

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
    
