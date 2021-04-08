include("objects.jl")

struct TaitEoS{T}
    B::T
    P₀::T
    C::T
    ρ₀::T
end

struct IdealGasEoS{T}
    T::T
    M::T
end

Broadcast.broadcastable(eos::TaitEoS) = Ref(eos)
Broadcast.BroadcastStyle(::Type{<:TaitEoS}) = 
                        Broadcast.DefaultArrayStyle{0}()

Broadcast.broadcastable(eos::IdealGasEoS) = Ref(eos)
Broadcast.BroadcastStyle(::Type{<:IdealGasEoS}) = 
                        Broadcast.DefaultArrayStyle{0}()

function density(x::TaitEoS, P::AbstractFloat)
    return x.ρ₀ / (1 - x.C * log10((x.B + P) 
                                   / (x.B + x.P₀)))
end

function density(x::IdealGasEoS, P::AbstractFloat)
    R = 8.314
    return P * x.M / R / x.T
end

function density!(sys::System)
    @. @views sys.ρ[2, :, :] = density(tait_C₅H₁₂, sys.P)
    @. @views sys.ρ[1, :, :] = density(ideal_gas, sys.P)
end

tait_C₅H₁₂ = TaitEoS(35e6, 1e5, 0.2105, 616.18)

ideal_gas = IdealGasEoS(298.0, 0.029)

function binarySearch(; f, left, right, niter=50, eps_x=1e-6)
    mid = left
    for i = 1 : niter
        mid = (left + right) / 2
        if abs(right - left) < eps_x
            return mid
        elseif f(mid) < 0
            left = mid
        elseif f(mid) > 0
            right = mid
        end
    end
    return error("root of function couldn't be found.")
end

function findPressure(; ρ̂₁::T, ρ̂₂::T, V::T,
        left=1e4, right=1e7, niter=50,
        eps_x=1e-6) where T<:AbstractFloat
    function f(P::Float64)
        ρ₁ = 0.029 / 8.314 / 298 * P
        s = ρ̂₁ / ρ₁
        ρ₂ = ρ̂₂ / (1 - s)
        ρ₀ = 616.18
        P₀ = 1e5
        return (ρ₂ - ρ₀) / ρ₂ - 0.2105 *
            log10((35e6 + P) / (35e6 + P₀))
    end

    P = binarySearch(; f, left, right, eps_x)
    ρ₁ = 0.029 / 8.314 / 298 * P
    s = ρ̂₁ / ρ₁
    return P, s
end

function findPressure!(sys::System, p::Parameters)
    for index in CartesianIndices(sys.P)
        i, j = Tuple(index)
        sys.P[i, j], sys.s[i, j] = findPressure(
                                    ρ̂₁=sys.ρ[1, i, j],
                                    ρ̂₂=sys.ρ[2, i, j],
                                    V=p.Δx * p.Δy)
    end
end
