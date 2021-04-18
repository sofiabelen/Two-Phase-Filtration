module SimulationParameters

export Δt, nsteps, L, nx, ny, Δx, Δy, φ, K, μ, Pin, Pout,
       P₀, M, ψ, ψ₀

# duration = 2000
Δt = 0.000001
# nsteps = round(Int64, duration / Δt)
nsteps = 0

## Space grid [0, 2] × [0, 2]
## Change later to 2 (how ?)
L = 2
nx = ny = 40
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
M = [0.028, 0.07215]

## Molar Composition on the Inlet
## ψ = ν₁ / ν₂ = (m₁ / M₁) / (m₂ / M₂)
ψ = 0.3

## Initial Molar Composition
ψ₀ = ψ

end
