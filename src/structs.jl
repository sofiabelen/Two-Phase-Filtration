## Index convention:
## [i, j, k], where k = 1, 2 is the component: 
## 1 - gas,
## 2 - liquid;
## [i, j] -> x, y

## -------------------------- BC -------------------------- ##
## (a f(x, y) + b ∂F / ∂n)|_{∂Ω} = g
## 
## We keep a and b in an unmutable struct because they are
## set as a parameter in the beginning of the simulation,
## and they don't change.
## Therefore, they are kept with the other simulation parameters
## 
## g, on the other hand, will get updated on every step,
## so it seems logical to keep together with all the other
## system parameters that change over time.
## We will denote it as f_RHS, since it is then easy to infer
## that we are refering the right hand side of the BC equation.
## 
## The struct BoundaryCondtions holds two 3D arrays, a and b.
## 
## Each of this holds for each component 4 arrays,
## one for each border: left, right, top, bottom.
## 
## This is stored in a 3D array, where the last index
## k corresponds to the component, and the second index
## denotes the border.
## 
## h[:][1][k] -> left
## 
## h[:][2][k] -> right
## 
## h[:][3][k] -> top
## 
## h[:][4][k] -> bottom
## 
## Pressure is the same for both components, but for 
## consistency we will also use the same structure, but
## use only the first index k = 1.

struct BoundaryCondition{T}
    a::Array{T, 3}
    b::Array{T, 3}
end

## Here we keep the parameters a and b for each variable.
struct BoundaryConditions{T}
    ## We don't need these thanks to our staggered grid
    # u::BoundaryCondition{T}
    # v::BoundaryCondition{T}
    ρ::BoundaryCondition{T}
    s::BoundaryCondition{T}
    P::BoundaryCondition{T}
end

## f_RHS == g
mutable struct BoundaryRHS{T}
    # u::Array{T, 3}
    # v::Array{T, 3}
    ρ::Array{T, 3}
    s::Array{T, 3}
    P::Array{T, 3}
end

Base.copy(bc_RHS::BoundaryRHS) = BoundaryRHS(copy.((## bc_RHS.u,
                                                    ## bc_RHS.v,
                                                    bc_RHS.ρ,
                                                    bc_RHS.s,
                                                    bc_RHS.P,))... )

## -------------------------------------------------------- ##

mutable struct System{T<:AbstractFloat}
    u::Array{T, 3}
    v::Array{T, 3}
    ρ::Array{T, 3}
    F::Array{T, 3}
    s::Array{T, 3}
    P::Array{T, 2}
    ## Mass flux along x axis
    Φu::Array{T, 3}
    ## Mass flux along x axis
    Φv::Array{T, 3}
    bc_RHS::BoundaryRHS{T}
end

Base.copy(sys::System) = System(copy.((sys.u, sys.v, sys.ρ,
                                       sys.F, sys.s, sys.P,
                                       sys.Φu, sys.Φv,
                                       sys.bc_RHS))... )

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
    bc::BoundaryConditions{T}
    wall::UnitRange{Int}
    inlet::UnitRange{Int}
    relax_steps::Integer
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
