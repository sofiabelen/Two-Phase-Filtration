## index [i, j, k], where k=1,2 is the component: 
## 1 - gas,
## 2 - liquid;
## [i, j] -> x, y
#
mutable struct System{T<:AbstractFloat}
    u::Array{T, 3}
    v::Array{T, 3}
    ρ::Array{T, 3}
    s::Array{T, 3}
    P::Array{T, 2}
end

Base.copy(sys::System) = System(copy.((sys.u, sys.v, sys.ρ,
                                       sys.s, sys.P))... )

struct Parameters{T<:AbstractFloat}
    nsteps::Integer
    nx::Integer
    ny::Integer
    L::Integer
    Δx::T
    Δy::T
    Δt::T
    φ::T
    K::T
    Pin::T
    Pout::T
    P₀::T
    ψ::T
    ψ₀::T
    M::Vector{T}
    μ::Vector{T}
end

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
