% Steady state PLG calculation
% Calculation to determine the number of bound PLG species

function [X,Cs_PLG] = quasiSteadyPLGCalculation(X, I, C_PLG_0)

global N_CS N_AV U_0 k_a k_r

C_f = N_AV*U_0*1E-15*1E-6;
k_an = k_a*N_CS/C_f;

C_0 = C_PLG_0*U_0*1E-6*1E-15*N_AV/N_CS;  

    for it_2 = 1:N_CS
        
        N = X(it_2 + 3*N_CS) - X(it_2) - X(it_2 + 2*N_CS);
        
        b = -(N + C_0 - I(it_2) + k_r(2)/k_an(2));
        woof_PLG = (-b - sqrt(b^2 - 4*(C_0 - I(it_2))*N))/2 ;
        
        if (woof_PLG) < 0
           
            woof_PLG = 0;
        end
        
        X(it_2 + 1*N_CS) = woof_PLG;
    
    end
   %hello = sum(X(N_CS   + 1:2*N_CS))/(CP.U_0*1E-6*1E-15*OP.N_AV/N_CS);

    %Cs_PLG = (C_0 - sum(X(N_CS   + 1:2*N_CS))*0.01)/(U_0*1E-6*1E-15*N_AV/N_CS)
    
    Cs_PLG = C_PLG_0 - (sum(X(N_CS   + 1:2*N_CS)) + sum(I))/(U_0*1E-6*1E-15*N_AV);
    
end