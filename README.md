# Hybrid stochastic-continuous model of fibrinolysis #

This folder contains the code required to generate Chapter 7 as laid out in the PhD thesis. 
For information on how the model works, please see the doumentation.

## Update history

Last updated May 2017

## Overview

This repository contains the following folders:

* engineSourceFiles: Main Matlab files that are used in the simulation. They contain the function files needed to run the hyrbid model.
* continuousKinetics: Files that solver for clot lysis in a continuous system (fibre shrinking model)
* stochasticShrinkModel: Files for the the hybrid model that solve fibrinolysis as the shrinking of the fibrin fibre network
* stochasticTransModel: Files for the hybrid model that solve fibrinolysis as the transverse cutting of fibrin fibres
* plottingFiles: Files that automate plotting of results.  

## Getting started ##

### Prerequisites
These files were created with Matlab 2013a and 2017a. 

### Deployment
Models can be called from Matlab scripts that call the source files needed to run. 

For instance, to run the continous model, check the funHybridShrinkModel.m file for the main function file. 

```
function [timespan,plotCon] = funHyrbidShrinkModel(R_f0,rho_0,finalTime,nt,C_tPA0,C_PLG0,C_PLS0)
```

The return arguments in this case is the time, and a structure that contains the variables of interest.  This can then be plotted 
in what ever way is needed. 

For instance, the compareDifferentFibrinolysisPatterns.m calls seperate models to compare their results:

```
fprintf('Running Hybrid model... \n')
tic()
[t,ploti] = funHybridTransModel(R_f0,rho_0,U_0,finalTime,nt,[Cs_tPA, Cs_PLG, Cs_PLS],N_fibres,lysisLevel);
toc()

% Running continuous method
fprintf('Running Continuous model... \n')
tic()
[timespan,plotCon] = funHybridShrinkModel(R_f0,rho_0,U_0,finalTime,nt,[Cs_tPA, Cs_PLG, Cs_PLS],N_fibres);
toc()

%save(strcat(num2str(U_0),'_comparison'));

plot(t,ploti.extentLysis,'-k',timespan,plotCon.extentLysis,'-r')
hold on
```

Note that inputs have to specified regarding the initial structure of the fibrin matrix

```
R_f0 = 100; % Intial radius of the fibrin network
rho_0 = 3;  % Initial concentration of the fibrin density 
U_0 = 100;  % Volume of the control envelope
finalTime = 10*60; % Final lysis time
nt = 1000;   % Number of timesteps
Cs_tPA = 0.04; Cs_PLG = 2.2; Cs_PLS = 0; % Initial concentrations of free phase proteins
lysisLevel = 0.2; % Critical level of lysis for the hybrid transmodel (nunber of remianing binding sites until the fibre is considered cleaved)
N_fibres = 1; % Tying fibres into packages to ease computation time (i.e. 100 would treate 100 crossections as 1 cross-section)
```

To examine the physics, check the main source files in engineSoruceFiles folder. 
Case specific functions can be found in their respective folders. 


