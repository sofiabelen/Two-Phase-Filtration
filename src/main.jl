include("cfd.jl")
include("objects.jl")
include("thermo.jl")
include("plot.jl")

using DelimitedFiles

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
    # duration = 2000
    Δt = 0.000001
    # nsteps = round(Int64, duration / Δt)
    nsteps = 10
    
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

    p = Parameters(nsteps, nx, ny, L, Δx, Δy, Δt, 
                        φ, K, Pin, Pout, P₀, ψ, ψ₀, M, μ)

    sys = filtration!(p)
    # plot(sys, p)
    # writedlm("pressure.txt", sys.P, ' ')
    # writedlm("u1.txt", sys.u[1, :, :], ' ')
    # writedlm("u2.txt", sys.u[2, :, :], ' ')
    # writedlm("density1.txt", sys.ρ[1, :, :], ' ')
    # writedlm("density2.txt", sys.ρ[2, :, :], ' ')
end
