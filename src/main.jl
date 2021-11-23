include("cfd.jl")
include("structs.jl")
include("thermo.jl")
include("dump.jl")
include("plot.jl")
include("test.jl")
include("SimulationParameters.jl")

## https://stackoverflow.com/questions/37200025/how-to-import-custom-module-in-julia
using DelimitedFiles, ..SimulationParameters
using Printf
sp = SimulationParameters


function update!(syswork::System, sysnext::System, p::Parameters,
    t::Integer)
# -------------------- Predictor ------------------------- #
## ρ̃ᵢⱼnext = ρᵢⱼwork + F⋅Δt      
    # sysnext_pred = copy(sysnext)

    calculate_flux!(syswork, p)
    continuity_equation!(syswork, sysnext, p)
    find_pressure!(sysnext, p)
    apply_bc_pres_sat_den!(sysnext, p, t)
    darcy!(sysnext, p)
    # apply_bc_velocity!(sysnext, p)
# -------------------------------------------------------- #
 
# # -------------------- Corrector ------------------------- #
# ## ρᵢⱼnext = ρᵢⱼwork + (F + F̃) / 2 ⋅Δt      
#     calculate_flux!(sysnext_pred, p)
#     continuity_equation!(syswork, sysnext_pred, sysnext, p)
#     find_pressure!(sysnext, p)
#     apply_bc_pres_sat_den!(sysnext, p, t)
#     darcy!(sysnext, p)
#     # apply_bc_velocity!(sysnext, p)
# # -------------------------------------------------------- #
end

function init(p::Parameters)
    check_bc(p.bc)
    P  = fill(p.P₀, p.nx, p.ny)
    s  = zeros(p.nx, p.ny, 2)
    u  = zeros(p.nx, p.ny, 2)
    v  = zeros(p.nx, p.ny, 2)
    ρ  = zeros(p.nx, p.ny, 2)
    F  = zeros(p.nx, p.ny, 2)
    Φu = zeros(p.nx, p.ny, 2)
    Φv = zeros(p.nx, p.ny, 2)
    bc_RHS = initial_RHS(p)

    sys = System(u, v, ρ, F, s, P, Φu, Φv, bc_RHS)

    initial_conditions!(sys, p)

    return sys
end

function filtration(p::Parameters)
    syswork = init(p)
    sysnext = copy(syswork)

    for t = 1 : p.nsteps
        check(syswork, sysnext, p, t)
        println("Step i = ", t)
        update!(syswork, sysnext, p, t)
        # dump_mass(syswork, sysnext, p, t)
        dump_total_in_out_flux(syswork, p, t)
        # dump(syswork, p)
        sysnext, syswork = syswork, sysnext
    end
    return syswork
end

function main(p)
    sys = filtration(p)
    plot(sys, p)
    dump(sys, p)
end

## -------------------------- BC -------------------------- ##

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

## Problem!!
## https://discourse.julialang.org/t/writing-functions-for-types-defined-in-another-module/31895/4
## https://discourse.julialang.org/t/structs-in-modules/34317
println("TODO: fix problem")
println("typeof(sp.bc): ", typeof(sp.bc))
println("typeof(bc): ", typeof(bc))

p = Parameters(sp.nsteps, sp.nx, sp.ny, sp.L, sp.Δx,
               sp.Δy, sp.Δt, sp.φ, sp.K, sp.Pin,
               sp.Pout, sp.P₀, sp.ψ, sp.ψ₀, sp.M, sp.μ,
               bc, sp.wall, sp.inlet, sp.relax_steps)

main(p)
