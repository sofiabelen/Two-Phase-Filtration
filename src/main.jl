include("cfd.jl")
include("structs.jl")
include("thermo.jl")
include("dump.jl")
include("SimulationParameters.jl")

## https://stackoverflow.com/questions/37200025/how-to-import-custom-module-in-julia
using DelimitedFiles, .SimulationParameters
sp = SimulationParameters

p = Parameters(sp.nsteps, sp.nx, sp.ny, sp.L, sp.Δx,
               sp.Δy, sp.Δt, sp.φ, sp.K, sp.Pin,
               sp.Pout, sp.P₀, sp.ψ, sp.ψ₀, sp.M, sp.μ)

function check_pressure(sys::System)
    check = 0
    for index in CartesianIndices(sys.P)
        i, j = Tuple(index)
        if sys.P[i, j] < 0
            check = 1
        end
    end
    if check == 0
        println("Pressure check: all good!")
    else
        println(sys.P)
        println("Pressure check: NEGATIVE PRESSURE!")
    end
end

function update!(syswork::System, sysnext::System,
        p::Parameters)
# -------------------- Predictor ------------------------- #
## ρ̃ᵢⱼnext = ρᵢⱼwork + F⋅Δt      
    sysnext_pred = copy(sysnext)

    continuity_equation!(syswork, sysnext_pred, p)
    find_pressure!(sysnext_pred, p)
    boundary_pressure!(sysnext_pred, p)
    boundary_density!(sysnext_pred, p) 
    boundary_saturation!(sysnext_pred, p)
    darcy!(sysnext_pred, p)
    boundary_velocity!(sysnext_pred, p)
# -------------------- Corrector ------------------------- #
## ρᵢⱼnext = ρᵢⱼwork + (F + F̃) / 2 ⋅Δt      
    continuity_equation!(syswork, sysnext_pred, sysnext, p)
    find_pressure!(sysnext, p)
    boundary_pressure!(sysnext, p)
    boundary_density!(sysnext, p) 
    boundary_saturation!(sysnext, p)
    darcy!(sysnext, p)
    boundary_velocity!(sysnext, p)
end

function init(p::Parameters)
    P = fill(p.P₀, p.nx, p.ny)
    s = zeros(p.nx, p.ny, 2)
    u = zeros(p.nx, p.ny, 2)
    v = zeros(p.nx, p.ny, 2)
    ρ = zeros(p.nx, p.ny, 2)
    F = zeros(p.nx, p.ny, 2)

    sys = System(u, v, ρ, F, s, P)

    boundary_pressure!(sys, p)
    initial_saturation!(sys, p)
    density!(sys)
    boundary_velocity!(sys, p)

    return sys
end

function filtration(p::Parameters)
    syswork = init(p)
    sysnext = copy(syswork)

    for t = 1 : p.nsteps
        println("Step i = ", t)
        update!(syswork, sysnext, p)
        sysnext, syswork = syswork, sysnext
    end
    return syswork
end

function main(p)
    sys = filtration(p)
    plot(sys, p)
    dump(sys, p)
end

main(p)
