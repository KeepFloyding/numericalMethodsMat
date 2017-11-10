function a = propensityFunctionPar3(X,Cs_tPA, Cs_PLG, Cs_PLS,epsilon)

global check

%==================================================================
% INITIALISATION
%==================================================================


%==================================================================
% DEPENDANT VARIABLES
%==================================================================

global k_a k_r k_2 k_cat K_M N_AV N_CS U_0 


C_f = N_AV*U_0*1E-15*1E-6;
k_an = k_a*N_CS/C_f;
K_Mn = K_M*C_f*(1/N_CS);

% Number of free phase molecules

Nfs_tPA = zeros(1,N_CS);

for it_tPA = 1 :N_CS
    Nfs_tPA(it_tPA) = (Cs_tPA*U_0*1E-6*1E-15*N_AV/N_CS);
end

Nfs_PLG = Cs_PLG*U_0*1E-6*1E-15*N_AV/N_CS;
Nfs_PLS = Cs_PLS*U_0*1E-6*1E-15*N_AV/N_CS;


% Number of free binding sites

N_free = zeros(1,N_CS);

for it_free = 1:N_CS
    
    N_free(it_free) = X(it_free + 3*N_CS) - X(it_free) - ...
        X(it_free + N_CS) - X(it_free + 2*N_CS);
    
    if (N_free(it_free) < 0)
        
        N_free(it_free) = 0;
        
    end
    
end

    for it = 1:length(X)
        
        if X(it)<0
            X(it) = 0;
        end
        
        
    end

%==================================================================
% GENERATING PROPENSITY FUNCTION
%==================================================================

% Number of rxns: 9
% Order is: PLS Dis, PLS Ads, M-M, cat, L_PLS Dis, tPA Ads, tPA Dis, PLG
% Ads, PLG Dis

    PLS_array = X(1 + 2*N_CS:3*N_CS);
    PLG_array = X(1 + 1*N_CS:2*N_CS);
    tPA_array = X(1:N_CS);
    c = zeros(1, N_CS); L_PLS_array = X(1 + 4*N_CS:5*N_CS);

    k_cat_new = k_cat;
    k_r_PLS = k_r(3); k_a_PLS = k_an(3);
    k_r_tPA = k_r(1); k_a_tPA = k_an(1);
    check_new = check;
    k_2_new = k_2;
             
        for it_G = 1:N_CS
            % M-M reaction
            c(it_G) = k_2_new*PLG_array(it_G)*X(it_G)/(K_Mn*(1-epsilon) + PLG_array(it_G));
            if (check_new == 1)
                
                c(it_G) = 0;
                
            end     
        end

    %clear PLS_array PLG_array L_PLS_array
    a = zeros(1, length(X));
    
    a(1:N_CS) = k_r_PLS*PLS_array; % PLS desorption
    a(N_CS+1:2*N_CS) = k_a_PLS*(Nfs_PLS/epsilon)*N_free; % PLS adsorption
    a(2*N_CS+1:3*N_CS) = c; % M-M reaction
    a(3*N_CS+1:4*N_CS) = k_cat_new*PLS_array; % Cat reaction
    a(4*N_CS+1:5*N_CS) = k_r_PLS*L_PLS_array; % L_PLS dis
    a(5*N_CS+1:6*N_CS) = k_a_tPA*(Nfs_tPA(it_G)/epsilon).*N_free; % tPA adsorption
    a(6*N_CS+1:7*N_CS) = k_r_tPA*tPA_array; % tPA desorption
    
    % PLG adsorption/desorption
    %a(it_G + 7*N_CS) = 0;%k_an(2)*(Nfs_PLG/epsilon)*N_free(it_G);
    %a(it_G + 8*N_CS) = 0;%k_r(2)*X(it_G + 1*N_CS);




end