% ======================================================================
% STOCHIOMETRY MATRIX
% ======================================================================

v = zeros(N_rxns,N_species);

v(1,3) = -1; % PLS desorption
v(2,3) =  1; % PLS adsorption

v(3,3)   = 1;  % M-M rxn
v(3,2)   = -1;

v(4,4)   = -1; % Cat rxn
v(4,3)   = -1;
v(4,5)   =  1;

v(5,5)   = -1; % L_PLS dis

v(6,1) = 1;             % tPA adsorption
v(7,1) = -1;            % tPA desorption

v(8,2) =  1; % PLG adsorption
v(9,2) = -1; % PLG desorption

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
    
