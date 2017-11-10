% Gillespie's Stochastic Method; the main file to edit would be the
% propensity function file. However, there is also a user-defined zone
% should any variables need to be computed outside the while loop.
% Otherwise, don't touch the other parts.

function [X_nx,Cs_PLG,I] = gillespieSolver(X,v, delta_t,propensityFunction,Cs_tPA, Cs_PLG, Cs_PLS,epsilon,I)
global N_CS epsilon_0 N_species

%==================================================================
% INITIALISING ALL VARIABLES
%==================================================================

t_G = 0;

%==================================================================
% USER-DEFINED ZONE
%==================================================================
% NOTE: Insert here any functions or calculation that need to be calculated
% outside the while loop to save on computational time.

% PLG Rapid SS Calculation

[X,Cs_PLG] = quasiSteadyPLGCalculation(X, I, 2.2*epsilon_0);


    % START OF WHILE LOOP
    while t_G < delta_t

        
        %==================================================================
        % ASSIGN RANDOM NUMBER (MONTE-CARLO STEP)
        %==================================================================
        r_1 = rand();
        r_2 = rand();
        
        %==================================================================
        % CALCULATE VALUE OF PROPENSITY FUNCTION
        %==================================================================
        
        a = propensityFunction(X,Cs_tPA, Cs_PLG, Cs_PLS,epsilon);
        
        %==================================================================
        % COMPUTE TIME-STEP AND REACTION INDEX
        %==================================================================
        
        
        a_0 = sum(a);
        dt = 1/a_0*log(1/r_1);  % timestep
        in = find(cumsum(a) > (a_0*r_2), 1 ); % reaction index
        
        if isempty(in);
            
            break;
        end
        
        
    if in >= 2*N_CS + 1 && in <= 3*N_CS 
       
        I(in - 2*N_CS) = I(in - 2*N_CS) + 1;
        
    end
        
        %==================================================================
        % UPDATE SIMULATION TIME 
        %==================================================================
        
        t_G = t_G + dt;
        
        if (t_G > delta_t)
            break;
        end
        
        %==================================================================
        % UPDATE STATE VECTOR
        %==================================================================
        
        CS_ind = rem(in,N_CS);
                
        if CS_ind == 0
           CS_ind = N_CS; 
        end

        rxn_ind = (in - CS_ind)/N_CS + 1;
                
        for it_inner = 1:N_species
            
            X((it_inner-1)*N_CS + CS_ind) = X((it_inner-1)*N_CS + CS_ind) + v(rxn_ind,it_inner);
                       
        end
        
        %X = X + v(in,:);
        
    end
    %END OF WHILE LOOP
    
    %==================================================================
    % ASSIGNING OUTPUT
    %==================================================================
    
    X_nx = X;

end