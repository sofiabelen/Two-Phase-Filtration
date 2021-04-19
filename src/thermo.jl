include("structs.jl")
include("physics.jl")

Broadcast.broadcastable(eos::TaitEoS) = Ref(eos)
Broadcast.BroadcastStyle(::Type{<:TaitEoS}) = 
                        Broadcast.DefaultArrayStyle{0}()

Broadcast.broadcastable(eos::IdealGasEoS) = Ref(eos)
Broadcast.BroadcastStyle(::Type{<:IdealGasEoS}) = 
                        Broadcast.DefaultArrayStyle{0}()

# -------------------- Density --------------------------- #

## The tait and ideal gas equations helps us determine
## ρ̂ᵢ= mᵢ/ sᵢV , but we need ρᵢ = mᵢ / V because ρᵢ is
## used in continuity and Darcy's equations. The first
## two functions return ρ̂ᵢ, while the third one returns ρᵢ

function density(x::TaitEoS, P::AbstractFloat)
    return x.ρ₀ / (1 - x.C * log10((x.B + P) 
                                   / (x.B + x.P₀)))
end

function density(x::IdealGasEoS, P::AbstractFloat)
    return P * x.M / R / x.T
end

function density!(sys::System)
    @. @views sys.ρ[:, :, 1] = density(ideal_gas, sys.P) *
                                      sys.s
    @. @views sys.ρ[:, :, 2] = density(tait_C₅H₁₂, sys.P) *
                                      (1 - sys.s)
end
# -------------------------------------------------------- #


# -------------------- Pressure -------------------------- #
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

function findPressure(; ρ₁::T, ρ₂::T, V::T,
        tait::TaitEoS, igas::IdealGasEoS,
        p::Parameters,
        left=1e4, right=1e7, niter=50,
        eps_x=1e-6) where T<:AbstractFloat

    function f(P::AbstractFloat)
        ρ̂₁ = density(igas, P)
        s = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return (ρ̂₂ - tait.ρ₀) / ρ̂₂ - tait.C *
            log10((tait.B + P) / (tait.B + tait.P₀))
    end

    P = binarySearch(; f, left, right, eps_x)
    ρ̂₁ = density(igas, P)
    s = ρ₁ / ρ̂₁

    return P, s
end

function findPressure!(sys::System, p::Parameters)
    for index in CartesianIndices(sys.P)
        i, j = Tuple(index)
        sys.P[i, j], sys.s[i, j] = findPressure(p=p,
                                    ρ₁=sys.ρ[i, j, 1],
                                    ρ₂=sys.ρ[i, j, 2],
                                    V=p.Δx * p.Δy,
                                    tait=tait_C₅H₁₂,
                                    igas=ideal_gas)
    end
end
# -------------------------------------------------------- #
