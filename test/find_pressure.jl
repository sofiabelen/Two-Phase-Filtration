include("../src/thermo.jl")
include("../src/physics.jl")
include("../src/SimulationParameters.jl")

using PyPlot

function plot_function(f, x₀)
    fig = PyPlot.figure(figsize=(10, 10))
    x = collect(-1e5:1e4:1e6)
    y = f.(x)
    grid("on")
    # vlines(1e8, -1000, 10000)
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
        ρ̂₁ = density(igas, P)
        s = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return (ρ̂₂ - tait.ρ₀) / ρ̂₂ - tait.C *
            log10((tait.B + P) / (tait.B + tait.P₀))
    end

    ρ₁ = 1.125
    ρ₂ = 2.898

    P, s = find_pressure(; ρ₁=ρ₁, ρ₂=ρ₂, p=p,
                        igas=ideal_gas, tait=tait_C₅H₁₂)

    println("P = ", P)
    println("s = ", s)

    plot_function(f, P)
end
