module SimulationParameters

include("physics.jl")

export Δt, nsteps, L, nx, ny, Δx, Δy, φ, K, μ, Pin, Pout,
       P₀, M, ψ, ψ₀

## Keep in mind: dt ∼ dx² !!

# Δt = 0.01
Δt = 0.001
# Δt = 0.0001
# Δt = 0.0004
# nsteps = 10000
# nsteps = 50000
# nsteps = 25000
nsteps = 50000

## Space grid [0, 2] × [0, 2]
L = 2
# nx = ny = 20
nx = ny = 50
# nx = ny = 150
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

## ψ = m₁ / (m₂ M₁ / M₂ + m₁)
## Molar Composition on the Inlet
ψ = 0.9

## Initial Molar Composition
ψ₀ = 0.1

end
