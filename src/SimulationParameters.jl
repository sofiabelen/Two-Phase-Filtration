module SimulationParameters

include("physics.jl")

export Δt, nsteps, L, nx, ny, Δx, Δy, φ, K, μ, Pin, Pout,
       P₀, M, ψ, ψ₀

Δt = 0.001
nsteps = 100000
# nsteps = 2

## Space grid [0, 2] × [0, 2]
L = 2
nx = ny = 50
Δx = Δy = L / nx

## Porosity
φ = 0.7

## Specific permeability
K = 1e-12

## Dynamic viscocity
μ = [1.8e-5, 2.14e-4]

## Pressures at top and bottom, and initial pressure
Pin = 1e6
Pout = 1e5
P₀ = Pout

## Molar masses
M = [MOLAR_MASS_IGAS, MOLAR_MASS_PENTANE]

## Molar Composition on the Inlet
ψ = 0.2

## Initial Molar Composition
ψ₀ = ψ

end
