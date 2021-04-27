include("../src/thermo.jl")
include("../src/physics.jl")
include("../src/SimulationParameters.jl")

using PyPlot

using .SimulationParameters
sp = SimulationParameters

p = Parameters(sp.nsteps, sp.nx, sp.ny, sp.L, sp.Δx,
               sp.Δy, sp.Δt, sp.φ, sp.K, sp.Pin,
               sp.Pout, sp.P₀, sp.ψ, sp.ψ₀, sp.M, sp.μ)

function plot_function(f::Function)
    fig = PyPlot.figure(figsize=(10, 10))
    x = collect(-1e5:1e4:1e6)
    y = f.(x)
    grid("on")
    # vlines(1e8, -1000, 10000)
    plot(x, y)
    show()
end

function plot_function(f::Function, g::Function)
    fig = PyPlot.figure(figsize=(10, 10))
    x = collect(-1e5:1e4:1e6)
    grid("on")

    plot(x, f.(x))
    plot(x, f.(x))
    plot(x, g.(x))

    show()
end

function test_find_pressure() 

    function f(P::AbstractFloat)
        igas = ideal_gas
        tait = tait_C₅H₁₂
        ρ̂₁ = density(igas, P)
        s = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return (ρ̂₂ - tait.ρ₀) / ρ̂₂ - tait.C *
            log10((tait.B + P) / (tait.B + tait.P₀))
    end

    function fder(P::AbstractFloat)
        igas=ideal_gas
        tait=tait_C₅H₁₂
        ρ̂₁ = density(igas, P)
        s = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return -tait.ρ₀ / ρ̂₂^2  * R * igas.T * igas.M *
            ρ₁ * ρ₂ / (igas.M * P - R * igas.T * ρ₁)^2 -
            tait.C / (tait.B + P) / log(10)
    end

    ρ₁ = 1.125
    ρ₂ = 2.898
    P, s = 0, 0

# -------------------- Binary Search --------------------- #
# 
#     println("Binary Search")
# 
#     @time begin
#     P, s = find_pressure(; root_finder=binary_search,
#                          ρ₁=ρ₁, ρ₂=ρ₂, p=p,
#                          igas=ideal_gas, tait=tait_C₅H₁₂)
#     end
# 
#     println("P = ", P)
#     println("s = ", s)
# -------------------------------------------------------- #

# -------------------- Secant Method --------------------- #
# 
#     println("Secant Method")
# 
#     @time begin
#     P, s = find_pressure(; root_finder=secant_root_finder,
#                          ρ₁=ρ₁, ρ₂=ρ₂, p=p,
#                          igas=ideal_gas, tait=tait_C₅H₁₂)
#     end
# 
#     println("P = ", P)
#     println("s = ", s)
# -------------------------------------------------------- #

    # plot_function(f, fder)
    # plot_function(fder)
# -------------------- Newton Raphson -------------------- #

    println("Newton Raphson")

    @time begin
    P, s = find_pressure(; root_finder=newton_raphson,
                         ρ₁=ρ₁, ρ₂=ρ₂, p=p,
                         igas=ideal_gas, tait=tait_C₅H₁₂)
    end

    println("P = ", P)
    println("s = ", s)
# -------------------------------------------------------- #


    # plot_function(f, P)
end

test_find_pressure()
