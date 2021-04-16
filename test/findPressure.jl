include("../src/thermo.jl")
include("../src/physics.jl")
include("../src/SimulationParameters.jl")

using PyPlot

function plotFunction(f, x₀)
    fig = PyPlot.figure(figsize=(10, 10))
    x = collect(500.0:1.0:1000.0)
    y = f.(x)
    grid("on")
    vlines(x₀, -.16, .125)
    plot(x, y)
    show()
end

let 
    using .SimulationParameters
    sp = SimulationParameters

    p = Parameters(sp.nsteps, sp.nx, sp.ny, sp.L, sp.Δx,
                   sp.Δy, sp.Δt, sp.φ, sp.K, sp.Pin,
                   sp.Pout, sp.P₀, sp.ψ, sp.ψ₀, sp.M, sp.μ)

    function f(P::AbstractFloat)
        igas = ideal_gas
        tait = tait_C₅H₁₂
        ρ₁ = density(igas, P)
        s = ρ̂₁ / ρ₁
        ρ₂ = ρ̂₂ / (1 - s)
        return (ρ₂ - tait.ρ₀) / ρ₂ - tait.C *
            log10((tait.B + P) / (tait.B + tait.P₀))
    end

    ρ̂₁ = 0.1
    ρ̂₂ = 1.0
    V = p.Δx * p.Δy

    P, s = (; ρ̂₁=ρ̂₁, ρ̂₂=ρ̂₂, V=V, p=p)

    println(P)
    println(s)

    plotFunction(f, P)
end
