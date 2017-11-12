% Plotting m-file

% =======================================================================
% INPUTS
% =======================================================================

lysisPercent = 0.05;

% =======================================================================
% FINDING LYSIS PARAMETERS
% =======================================================================

% lysis index and lysis line
lysisIndex = find(n_tot/n_tot(1) <= lysisPercent,1);
tLysis = t(lysisIndex);

% Line of lysis
lysisLine = ones(1,nt)*lysisPercent*n_tot(1);

% Converting to minutes
t_min = t/60; % Plotting in minutes
tLysisMin = tLysis/60;

figure(1)

% =======================================================================
% BOUND PHASE PROTEINS
% =======================================================================

subplot(2,2,1)
plot(t_min,n_tPA,t_min,n_PLG,t_min,n_PLS)

% =======================================================================
% VOIDAGE
% =======================================================================

subplot(2,2,2)
plot(t_min, epsilon,tLysisMin,epsilon(lysisIndex),'s')

% =======================================================================
% TOTAL BINDING SITES
% =======================================================================

subplot(2,2,3)
plot(t_min, n_tot,t_min, lysisLine)
hold on

% =======================================================================
% FREE PHASE PROTEINS
% =======================================================================

subplot(2,2,4)
plot(t_min,C_tPA,t_min,C_PLS,tLysisMin,C_PLS(lysisIndex),'s')
