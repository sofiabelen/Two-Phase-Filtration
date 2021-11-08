module SimulationParameters

include("physics.jl")

export Δt, nsteps, L, nx, ny, Δx, Δy, φ, K, μ, Pin, Pout,
       P₀, M, ψ, ψ₀, bc, inlet, wall

## CFL condition: dt ∼ dx² !!

Δt = 0.01
# Δt = 0.005
# Δt = 0.0001
# Δt = 0.0001
# nsteps = 10000
# nsteps = 50000
# nsteps = 25000
nsteps = 1000

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

## ψ = m₁ / (m₂ M₁ / M₂ + m₁)
## Molar Composition on the Inlet
ψ = 0.9

## Initial Molar Composition
ψ₀ = 0.1

## -------------------------- BC -------------------------- ##

inlet = 1 : nx ÷ 2
wall  = nx ÷ 2 + 1: nx

a_u = zeros(max(nx, ny), 4, 2)
a_v = zeros(max(nx, ny), 4, 2)
a_ρ = zeros(max(nx, ny), 4, 2)
a_s = zeros(max(nx, ny), 4, 2)
a_P = zeros(max(nx, ny), 4, 2)
                        
b_u = zeros(max(nx, ny), 4, 2)
b_v = zeros(max(nx, ny), 4, 2)
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

## top & bottom
@views a_v[:, :, k]     .= 1

## left & right
@views a_u[:, :, k]     .= 1

## left & right & top & bottom
@views a_ρ[:, :, k]     .= 1

## inlet at bottom [1, nx ÷ 2]
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

u_bc = BoundaryCondition(a_u, b_u)
v_bc = BoundaryCondition(a_v, b_v)
ρ_bc = BoundaryCondition(a_ρ, b_ρ)
s_bc = BoundaryCondition(a_s, b_s)
P_bc = BoundaryCondition(a_P, b_P)

function check_bc(bc::BoundaryCondition, var::String)
    for k = 1 : 2
        for w = 1 : 4
            for x = 1 : min(nx, ny)
                if (bc.a[x, w, k]^2 + bc.b[x, w, k]^2) == 0
                    return error("Wrong ", var, " bc at k = ", k,
                                 " w = ", w, " x = ", x)
                end
            end
        end
    end
end

check_bc(u_bc, "u")
check_bc(v_bc, "v")
check_bc(ρ_bc, "ρ")
check_bc(s_bc, "s")
check_bc(P_bc, "P")

bc = BoundaryConditions(u_bc, v_bc, ρ_bc, s_bc, P_bc)
end
