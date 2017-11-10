# Numerical Methods in Matlab

This repository contains Matlab functions of several popular numerical methods for solving differential equations such as Gillespie's stochastic method and Runge-Kutta 4th order method. 

## Getting started

### Prerequisites

* Matlab 2013a or higher (tested until Matlab 2017a)

### Installing 

Simply copy the Matlab files in src to your working directory, make sure that it is included in your path and call as required. 

### Testing 

A few test cases have been set up in test that showcases the use of the models in solving differential equations. Simply run runAll.m in Matlab console. 

## Deploying

### Gillespie stochastic method

You will need to edit the following files

* stochiometryMatrix.m as per the stochiometry of your reactions
* propensityFunction.m as per the reaction rate of each reaction and species. 

Specify an intial state array (X), a timestep (dt), stochiometry matrix (v) and a propensity function. 

Make sure to update the state array after every timestep.

```
X = zeros(1,N_species,'single')

for it = 2:nt
	X_nx = gillespieSolver(X,v, delta_t,@propensityFunction);
	X = X_nx;
end
```

### Runge-Kutta method

You will need to edit

* gradient.m as per the differential equations that you wish to solve. 

For instance

```
eqns(1)
```

represents the derivative of the first component with respect to time. 

In your main Matlab file, specify a timestep (delta_t) and the initial value of your variables (@ time = 0). Put the main equation in a for loop and update the time after ever timestep. 

```
input = zeros(1,N_eqns)

for it = 2:nt
    
    input = [x; y; z];
    current_t = (it-1)*delta_t;
    
    sol = rungeKuttaFourth(input, current_t, delta_t);

 end
```


