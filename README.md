# Hydrodynamics Solver: Ideal Gas and Pentane Filtration

## How to Run

In normal mode:

```
julia src/main.jl
```

In debug mode:

```
JULIA_DEBUG=Main julia src/main.jl
```

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

## 2-Component Filtration: Ideal Gas and Pentane
![Alt text](img/2phase-filtration.svg?raw=true "1")

## 1-Component Cavity Flow: Ideal Gas
![Alt text](img/ideal-gas.svg?raw=true "1")

## 1-Component Filtration: Ideal Gas
![Alt text](img/filtration.svg?raw=true "1")
