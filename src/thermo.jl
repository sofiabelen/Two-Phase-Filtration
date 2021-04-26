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
    if x.B + P < 0
        println("Less than 0: ", P)
    end
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
function binary_search(; f, left, right, niter=50, eps_x=1e-6)
    fleft = f(left)
    fright = f(right)
    mid = (left + right) / 2
    for i = 1 : niter
        mid = (left + right) / 2
        if abs(right - left) < eps_x
            return mid
        end
        fmid = f(mid)
        fmid == 0 && return mid
        if sign(fmid) == sign(fleft)
            left = mid
            fleft = fmid
        else
            right = mid
            fright = fmid
        end
    end
    return error("root of function couldn't be found: mid = ", mid)
end

function find_pressure(; ρ₁::T, ρ₂::T,
        tait::TaitEoS, igas::IdealGasEoS,
        p::Parameters, left=1e4, right=1e8, niter=50,
        eps_x=1e-6) where T<:AbstractFloat

    function f(P::AbstractFloat)
        ρ̂₁ = density(igas, P)
        s = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return (ρ̂₂ - tait.ρ₀) / ρ̂₂ - tait.C *
            log10((tait.B + P) / (tait.B + tait.P₀))
    end

    P = binary_search(; f, left, right, eps_x)
    ρ̂₁ = density(igas, P)
    s = ρ₁ / ρ̂₁

    return P, s
end

function find_pressure!(sys::System, p::Parameters)
    for index in CartesianIndices(sys.P)
        i, j = Tuple(index)
        sys.P[i, j], sys.s[i, j] = find_pressure(p=p,
                                    ρ₁=sys.ρ[i, j, 1],
                                    ρ₂=sys.ρ[i, j, 2],
                                    tait=tait_C₅H₁₂,
                                    igas=ideal_gas)
    end
    @debug "Pressure " sys.P
end
# -------------------------------------------------------- #
