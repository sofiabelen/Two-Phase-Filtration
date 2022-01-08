# Hydrodynamics Solver: Nitrogen and Pentane Filtration

Check out this [post](https://sofiabelen.github.io/projects/two-phase-filtration.html) for a more detailed description!

## Directory Structure

- `doc/` - compact documentation and progress tracking.
- `src/` - source code.
    - `main` - entry point.
    - `physics` - physical and EoS parameters.
    - `cfd` - hydrodynamics solver.
    - `thermo` - thermodynamics functions.
    - `computational` - computational methods.
    - `dump` - data dump + plot.
    - `SimulationParameters` - simulation parameters.
- `test/` - for unit testing.
- `aux/` - supplementary problems leading up to main project.

## Initial and Boundary Conditions
![Alt text](img/problem-formulation.svg?raw=true "1")

## Methods Used

![Alt text](img/staggered-grid.png?raw=true "1")

1. Second order finite difference method for spatial discretization using a staggered grid.

2. Explicit predictor-corrector method according to the Heun scheme  for time integration.

3. Newton-Raphson method for finding pressure and gas saturation.

## Algorithm

1. We are given the densities $\rho_i = \rho_i(t)$.

2. Finding the pressure and gas saturation using the
Newton-Raphson method from the condition of equality
of the pressure of the gas and the liquid:

$$P = P_1 \left( \frac{\rho_1}{s} \right) 
= P_2 \left( \frac{\rho_2}{1 - s}\right)$$

3. Calculation of fluxes from Darcy's law.

4. Calculation of the densities $\rho_i(t + \Delta t)$ based on the known fluxes.

5. Renaming $\rho_i = \rho_i(t + \Delta t)$ and moving on to the next time step.

## Results

![Alt text](img/2phase-filtration-density.png?raw=true "1")

![Alt text](img/fluxes.gif?raw=true "1")

## How to Run

TODO: PkgTemplates.jl

In normal mode:

```
julia src/main.jl
```

In debug mode:

```
JULIA_DEBUG=Main julia src/main.jl
```
