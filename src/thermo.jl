include("structs.jl")
include("physics.jl")
include("computational.jl")

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

function find_pressure(; root_finder::Function, ρ₁::T, ρ₂::T,
        tait::TaitEoS, igas::IdealGasEoS,
        p::Parameters, niter::Integer=50,
        eps_x::T=1e-6)::Tuple{T, T} where T<:AbstractFloat

    function f(P::T)::T where T<:AbstractFloat
        ρ̂₁ = density(igas, P)
        s  = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return (ρ̂₂ - tait.ρ₀) / ρ̂₂ - tait.C *
            log10((tait.B + P) / (tait.B + tait.P₀))
    end

    function fder(P::T)::T where T<:AbstractFloat
        ρ̂₁ = density(igas, P)
        s  = ρ₁ / ρ̂₁
        ρ̂₂ = ρ₂ / (1 - s)
        return -tait.ρ₀ / ρ̂₂^2  * R * igas.T * igas.M *
            ρ₁ * ρ₂ / (igas.M * P - R * igas.T * ρ₁)^2 -
            tait.C / (tait.B + P) / log(10)
    end

    let
        P  = root_finder(; f=f, fder=fder, eps_x=eps_x)
        ρ̂₁ = density(igas, P)
        s  = ρ₁ / ρ̂₁

        return P, s
    end
end

function find_pressure!(sys::System, p::Parameters)
    for index in CartesianIndices(sys.P)
        i, j = Tuple(index)
        ## This gives seg fault, but is type stable
        # (sys.P[i, j], sys.s[i, j])::Tuple{T, T} =
        #             find_pressure(p=p,
        #             root_finder=newton_raphson,
        #             ρ₁=sys.ρ[i, j, 1],
        #             ρ₂=sys.ρ[i, j, 2],
        #             tait=tait_C₅H₁₂,
        #             igas=ideal_gas) where T<:AbstractFloat
        #

        ## Type unstable
        sys.P[i, j], sys.s[i, j] =
                   find_pressure(p=p,
                   root_finder=newton_raphson,
                   ρ₁=sys.ρ[i, j, 1],
                   ρ₂=sys.ρ[i, j, 2],
                   tait=tait_C₅H₁₂,
                   igas=ideal_gas)
    end
end
# -------------------------------------------------------- #
