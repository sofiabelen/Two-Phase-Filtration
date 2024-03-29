module SimulationParameters

include("physics.jl")

export Δt, nsteps, L, nx, ny, Δx, Δy, φ, K, μ, Pin, Pout,
       P₀, M, ψ, ψ₀, bc, inlet, wall, relax_steps

## CFL condition: dt ∼ dx² !!

Δt = 0.01
# Δt = 0.005
# Δt = 0.0001
# Δt = 0.0001
# nsteps = 10000
# nsteps = 50000
# nsteps = 25000

# nsteps = 500000
# nsteps = 320000
nsteps = 1000

relax_t = 100
relax_steps = round(Int, relax_t / Δt)
println("relax steps: ", relax_steps)
# relax_steps = 100

## Space grid [0, 2] × [0, 2]
L = 2
nx = ny = 20
# nx = ny = 30
# nx = ny = 50
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

## ψ = ν₁ / ν₂ = m₁ / (m₂ M₁ / M₂ + m₁)
## Molar Composition on the Inlet
ψ = 0.9

## Initial Molar Composition
ψ₀ = 0.1

## -------------------------- BC -------------------------- ##

boundary = 1 : nx

inlet = 1 : nx ÷ 2
wall  = nx ÷ 2 + 1: nx
# inlet = nx ÷ 3 + 1: 2 * nx ÷ 3
# inlet = 1 : nx
# wall  = nx: nx
# This works as if we had an array of 'inlet'
# the same size of boundary
#wall = boundary[.!in.(boundary, Ref(inlet))]

## mask : (1:10)[rand(Bool, 10)]

# a_u = zeros(max(nx, ny), 4, 2)
# a_v = zeros(max(nx, ny), 4, 2)
a_ρ = zeros(max(nx, ny), 4, 2)
a_s = zeros(max(nx, ny), 4, 2)
a_P = zeros(max(nx, ny), 4, 2)
                        
# b_u = zeros(max(nx, ny), 4, 2)
# b_v = zeros(max(nx, ny), 4, 2)
b_ρ = zeros(max(nx, ny), 4, 2)
b_s = zeros(max(nx, ny), 4, 2)
b_P = zeros(max(nx, ny), 4, 2)

## h[:][1][k] -> left
## 
## h[:][2][k] -> right
## 
## h[:][3][k] -> top
## 
## h[:][4][k] -> bottom

k = 1 : 2

## velocities all
# @views a_v[:, :, k]     .= 1
# @views a_u[:, :, k]     .= 1

## left & right & top & bottom
@views a_ρ[:, :, k]     .= 1

## inlet at bottom [1, nx ÷ 2]
for m = 1 : n_in
@views a_s[inlet, 4, k] .= 1

## walls & outlet
@views b_s[:, 1:3,  k]  .= 1
@views b_s[wall, 4, k]  .= 1

## outlet
@views a_P[:, 3, k]     .= 1

## inlet
@views a_P[inlet, 4, k] .= 1

## walls left & right
@views b_P[:, 1:2, k]   .= 1

## bottom wall
@views b_P[wall, 4, k]  .= 1

# u_bc = BoundaryCondition(a_u, b_u)
# v_bc = BoundaryCondition(a_v, b_v)
ρ_bc = BoundaryCondition(a_ρ, b_ρ)
s_bc = BoundaryCondition(a_s, b_s)
P_bc = BoundaryCondition(a_P, b_P)

bc = BoundaryConditions(ρ_bc, s_bc, P_bc)
end
