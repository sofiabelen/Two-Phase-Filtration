include("cfd.jl")
include("structs.jl")
include("thermo.jl")
include("plot.jl")
include("SimulationParameters.jl")

function update!(syswork::System, sysnext::System,
        p::Parameters)
    continuity_equation!(syswork, sysnext, p)
    find_pressure!(sysnext, p)
    boundary_pressure!(sysnext, p)
    boundary_density!(sysnext, p) 
    boundary_saturation!(sysnext, p)
    darcy!(sysnext, p)
    boundary_velocity!(sysnext, p)
end

function init(p::Parameters)
    P = fill(p.P₀, p.nx, p.ny)
    s = zeros(p.nx, p.ny)
    u = zeros(p.nx, p.ny, 2)
    v = zeros(p.nx, p.ny, 2)
    ρ = zeros(p.nx, p.ny, 2)

    sys = System(u, v, ρ, P, s)

    boundary_pressure!(sys, p)
    initial_saturation!(sys, p)
    density!(sys)
    boundary_velocity!(sys, p)

    return sys
end

function filtration!(p::Parameters)
    syswork = init(p)
    sysnext = copy(syswork)

    for t = 1 : p.nsteps
        update!(syswork, sysnext, p)
        sysnext, syswork = syswork, sysnext
    end
    return syswork
end

let
    ## https://stackoverflow.com/questions/37200025/how-to-import-custom-module-in-julia
    using DelimitedFiles, .SimulationParameters
    sp = SimulationParameters

    p = Parameters(sp.nsteps, sp.nx, sp.ny, sp.L, sp.Δx,
                   sp.Δy, sp.Δt, sp.φ, sp.K, sp.Pin,
                   sp.Pout, sp.P₀, sp.ψ, sp.ψ₀, sp.M, sp.μ)

    sys = filtration!(p)
    plot(sys, p)
    x = ideal_gas
    pressure1 = @. R * x.T * sys.ρ[:, :, 1] / p.M[1] / sys.s
    x = tait_C₅H₁₂
    ρ̂₂ = @. sys.ρ[:, :, 2] / (1 - sys.s)
    pressure2 = @. (x.B + x.P₀) * 10^((ρ̂₂ - x.ρ₀) / x.C / 
                                      ρ̂₂) - x.B
    writedlm("pressure.txt", sys.P, ' ')
    writedlm("pressure1.txt", pressure1 .- sys.P, ' ')
    writedlm("pressure2.txt", pressure2 .- sys.P, ' ')
    writedlm("v1.txt", sys.v[:, :, 1], ' ')
    writedlm("v2.txt", sys.v[:, :, 2], ' ')
    writedlm("density1.txt", sys.ρ[:, :, 1], ' ')
    writedlm("density2.txt", sys.ρ[:, :, 2], ' ')
    writedlm("saturation.txt", sys.s, ' ')
end
