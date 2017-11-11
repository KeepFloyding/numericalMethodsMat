# Numerical Methods in Matlab

This repository contains Matlab functions of several popular numerical methods for solving differential equations such as Gillespie's stochastic method and Runge-Kutta 4th order method. 

## Getting started

### Prerequisites

* Matlab 2013a or higher (tested until Matlab 2017a)

### Installing 

Simply copy the Matlab files in src to your working directory, make sure that the file is included in your path and call as required. 

### Testing 

A few test cases have been set up in test that showcases the use of the models in solving differential equations. For instance, to test Runge-Kutta, simply run runTestRK.m in Matlab console. 

runTestAll.py runs both the Runge-Kutta test file and the Gillespie test file and compares the results with a plot. Both solvers solve the concentrations of species A, B and C that undergo the reversible reaction

```
A + B <-> C
```

You should get something like the following image where the smooth lines are the results calculated by Runge-Kutta methods and the erratic lines are those solved by the Gillespie method.

![testimage](https://user-images.githubusercontent.com/29730122/32693828-2d8bce30-c729-11e7-9dbe-b9432139ab9c.png)


## Deploying

### Gillespie stochastic method

This method requires the following:

* Intial state array (X)
* a timestep (dt)
* stochiometry matrix (v)
* propensityFunction.m as per the reaction rate of each reaction and species. 

In reality the model can simply be called as 

```
X_nx = gillespieSolver(X,v, final_time, @propensityFunction);
```

to get the value of the state array at the specified final_time. Since Gillespie's method gives a non-uniform time according to the values of the propensity function, it can be editted as follows to save values at every timestep dt. 

```
X = [N_A, N_B, N_C]

for it = 2:nt
	X_nx = gillespieSolver(X,v, delta_t,@propensityFunction);
	X = X_nx;
end
```

Make sure to update the state array after every timestep.


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
input = [x; y; z];

for it = 2:nt
    
    current_t = (it-1)*delta_t;
    
    sol = rungeKuttaFourth(input, current_t, delta_t);

    % Update the input
    input = [sol(1);sol(2);sol(3)]

 end
```


