using PyPlot

function binarySearch(; f, left, right, niter=50, eps=1e-6)
    mid = left
    for i = 1 : niter
        mid = left + (right - left) / 2

        if abs((right - left) / mid) < eps
            return mid
        elseif f(mid) < 0
            left = mid
        elseif f(mid) > 0
            right = mid
        else
            return mid
        end
    end
    return mid
end

function plotFunction(f, x₀)
    fig = PyPlot.figure(figsize=(10, 10))
    x = collect(500:1:1000)
    y = f.(x)
    grid("on")
    vlines(x₀, -.16, .125)
    plot(x, y)
    show()
end

function findPressure(; left, right, niter, m₁, m₂, V, eps)
    ## We want to find "any" zero of this function
    function f(ρ₂)
        ρ₁ = m₁ * ρ₂ / (V * ρ₂ - m₂)
        P = 8.314 * 298 / 0.029 * ρ₁
        ρ₀ = 616.18
        P₀ = 1e5
        return (ρ₂ - ρ₀) / ρ₂ - 0.2105 * log10((35e6 + P) / (35e6 + P₀))
    end

    x₀ = binarySearch(; f, left, right, eps)
    plotFunction(f, x₀)

    ρ₁ = m₁ * x₀ / (V * x₀ - m₂)
    println("ρ₁", ρ₁)
    P = 8.314 * 298 / 0.029 * ρ₁
    return P
end

let
    args = (
            left = 500,
            right = 1000,
            m₁ = 0.02,
            m₂ = 0.1,
            V = 0.01
            eps = 1e-6
           )
    println(findPressure(; args...))
end
