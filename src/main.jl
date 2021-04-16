include("cfd.jl")
include("structs.jl")
include("thermo.jl")
include("plot.jl")
include("SimulationParameters.jl")

function update!(syswork::System, sysnext::System,
        p::Parameters)
    continuity_equation!(syswork, sysnext, p)
    findPressure!(sysnext, p)
    boundaryPressure!(sysnext, p)
    boundaryDensity!(sysnext, p) 
    boundarySaturation!(sysnext, p)
    darcy!(sysnext, p)
    boundaryVelocity!(sysnext, p)
end

function init(p::Parameters)
    P = fill(p.P₀, p.nx, p.ny)
    u = zeros(2, p.nx, p.ny)
    v = zeros(2, p.nx, p.ny)
    s = zeros(p.nx, p.ny)
    ρ = zeros(2, p.nx, p.ny)

    sys = System(u, v, ρ, P, s)

    boundaryPressure!(sys, p)
    density!(sys)
    initialSaturation!(sys, p)
    boundaryVelocity!(sys, p)

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
    # plot(sys, p)
    # writedlm("pressure.txt", sys.P, ' ')
    # writedlm("u1.txt", sys.u[1, :, :], ' ')
    # writedlm("u2.txt", sys.u[2, :, :], ' ')
    # writedlm("density1.txt", sys.ρ[1, :, :], ' ')
    # writedlm("density2.txt", sys.ρ[2, :, :], ' ')
end
