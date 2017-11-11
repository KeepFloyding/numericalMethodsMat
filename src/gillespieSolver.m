% Gillespie's Stochastic Method; the main file to edit would be the
% propensity function file. However, there is also a user-defined zone
% should any variables need to be computed outside the while loop.
% Otherwise, don't touch the other parts.

function [X_nx,Cs_PLG,I] = gillespieSolver(X,v, delta_t,propensityFunction)


%==================================================================
% INITIALISING ALL VARIABLES
%==================================================================

t_G = 0;

%==================================================================
% SIMULATION ENGINE
%==================================================================


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
        
        a = propensityFunction(X);
        
        %==================================================================
        % COMPUTE TIME-STEP AND REACTION INDEX
        %==================================================================
        
        
        a_0 = sum(a);
        dt = 1/a_0*log(1/r_1);  % timestep
        in = find(cumsum(a) > (a_0*r_2), 1 ); % reaction index
        
        if isempty(in);
            
            break;
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
        
        X = X + v(in,:);
        
    end
    %END OF WHILE LOOP
    
    %==================================================================
    % ASSIGNING OUTPUT
    %==================================================================
    
    X_nx = X;

end