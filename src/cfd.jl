include("structs.jl")
include("thermo.jl")

## Avoid typos when working with partial derivatives

## NW--N--NE
##     |
## W---C---E
##     |
## SW--S--SE

function cnswe(arr::AbstractMatrix{T},
        i::Int, j::Int) where T<:AbstractFloat
    c = arr[i, j]
    n = arr[i, j + 1]
    s = arr[i, j - 1]
    w = arr[i - 1, j]
    e = arr[i + 1, j]
    return c, n, s, w, e
end

# ------------------ Continuity Equation ----------------- #

## Continuity equation to calculate ρᵢ
## φ ⋅∂ρᵢ / ∂t + div(ρᵢv⃗ᵢ) = 0, where ρᵢ = mᵢ / V
## We will add a right side to this where we inject matter.
function continuity_equation!(;
        ρnext::AbstractMatrix{T},
        ρwork::AbstractMatrix{T},
        uwork::AbstractMatrix{T},
        vwork::AbstractMatrix{T},
        p::Parameters) where T<:AbstractFloat
    for j = 2 : p.ny - 1
        for i = 2 : p.nx - 1
            ρc, ρn, ρs, ρw, ρe = cnswe(ρwork, i, j)
            uc, un, us, uw, ue = cnswe(uwork, i, j)
            vc, vn, vs, vw, ve = cnswe(vwork, i, j)
    
            ρnext[i, j] = ρc - p.Δt / p.φ *
                (uc * (ρe - ρw) / (2 * p.Δx) +
                 vc * (ρn - ρs) / (2 * p.Δy) +
                 ρc * (ue - uw) / (2 * p.Δx) +
                 ρc * (vn - vs) / (2 * p.Δy))
        end
    end
end

function continuity_equation!(syswork::System, 
        sysnext::System, p::Parameters)
    for k = 1 : 2
        @views continuity_equation!(ρnext=sysnext.ρ[:, :, k], 
                                    ρwork=syswork.ρ[:, :, k],
                                    uwork=syswork.u[:, :, k],
                                    vwork=syswork.v[:, :, k],
                                    p=p)
    end
end
# -------------------------------------------------------- #


# ---------------------- Darcy --------------------------- #

## Darcy's law
## v⃗ = -μ⁻¹ K̂⋅f(s) ⋅∇P

function darcy!(; f::Function,
        u::AbstractMatrix{T}, v::AbstractMatrix{T},
        P::AbstractMatrix{T}, s::AbstractMatrix{T},
        p::Parameters, k::Integer) where T<:AbstractFloat
    for j = 2 : p.ny - 1
        for i = 2 : p.nx - 1
            ## можно через (ρwork + ρ) / 2 -> 2-го порядка
            Pc, Pn, Ps, Pw, Pe = cnswe(P, i, j)

            u[i, j] = -p.K / p.μ[k] * f(s[i, j]) *
                (Pe - Pw) / (2 * p.Δx)
                v[i, j] = -p.K / p.μ[k] * f(s[i, j]) *
                (Pn - Ps) / (2 * p.Δy)
        end
    end
end

function darcy!(sys::System, p::Parameters)
    ## s is the saturation of gas
    fgas(s) = s^2
    fliq(s) = (1 .- s)^2
    f = [fgas, fliq]

    for k = 1 : 2
        @views darcy!(f=f[k], u=sys.u[:, :, k],
                      v=sys.v[:, :, k], P=sys.P,
                      s=sys.s, p=p, k=k)
    end
end
# -------------------------------------------------------- #



# ---------------------- BC & IC ------------------------- #

function boundary_density!(ρ::AbstractArray{T, 3},
        P::AbstractMatrix{T}, s::AbstractMatrix{T},
        p::Parameters) where T<:AbstractFloat
    @. @views ρ[:, 1, 1]  = density(ideal_gas, P[:, 1]) *
                                               s[:, 1]
    @. @views ρ[:, p.ny, 1] = density(ideal_gas, P[:, p.ny])*
                                                 s[:, p.ny]
    @. @views ρ[1, :, 1]  = density(ideal_gas, P[1, :]) *
                                               s[1, :]
    @. @views ρ[p.nx, :, 1] = density(ideal_gas, P[p.nx, :])*
                                                 s[p.nx, :]  

    @. @views ρ[:, 1, 2]  = density(tait_C₅H₁₂, P[:, 1]) *
                                           (1 - s[:, 1])
    @. @views ρ[:, p.ny, 2] = density(tait_C₅H₁₂, P[:,p.ny])*
                                             (1 - s[:,p.ny])
    @. @views ρ[1, :, 2]  = density(tait_C₅H₁₂, P[1, :]) *
                                           (1 - s[1, :])
    @. @views ρ[p.nx, :, 2] = density(tait_C₅H₁₂, P[p.nx,:])*
                                             (1 - s[p.nx,:])
end

function boundary_density!(sys::System, p::Parameters)
    boundary_density!(sys.ρ, sys.P, sys.s, p)
end

## Inlet
## We derive the saturation from the molar composition
## ψ = ν₁ / ν₂ and the densities ρ̂₁ and ρ̂₂.

## Outlet
## ∂s / ∂n⃗ = 0
function boundary_saturation!(
        ρ̂₁::AbstractMatrix{T},
        ρ̂₂::AbstractMatrix{T},
        s::AbstractMatrix{T},
        p::Parameters) where T<:AbstractFloat

    @. @views s[:, 1] = (ρ̂₁[:, 1] / ρ̂₂[:, 1] *
                         p.M[2] * (1 - p.ψ) / p.M[1] /
                         p.ψ + 1)^(-1)
    @. @views s[:, ny] = s[:, ny - 1]
end

function boundary_saturation!(sys::System, p::Parameters)
    ρ̂₁ = density.(ideal_gas, sys.P)
    ρ̂₂ = density.(tait_C₅H₁₂, sys.P)
    boundary_saturation!(ρ̂₁, ρ̂₂, sys.s, p)
end

function initial_saturation!(
        ρ̂₁::AbstractMatrix{T},
        ρ̂₂::AbstractMatrix{T},
        s::AbstractMatrix{T}, p) where T<:AbstractFloat

    @. @views s[:, :] = (ρ̂₁[:, :] /
                         ρ̂₂[:, :] * p.M[2] * (1 - p.ψ₀) /
                         p.M[1] / p.ψ₀ + 1)^(-1)
end

function initial_saturation!(sys::System, p::Parameters)
    ρ̂₁ = density.(ideal_gas, sys.P)
    ρ̂₂ = density.(tait_C₅H₁₂, sys.P)
    initial_saturation!(ρ̂₁, ρ̂₂, sys.s, p)
end

function boundary_pressure!(P::Array{T, 2},
        p::Parameters) where T<:AbstractFloat
    ## P = Pin at y = 0 and x ∈ [0, 1]
    @. @views P[:, 1] = 2 * p.Pin - P[:, 2]

    ## P = Pout at y = 2
    @. @views P[:, p.ny] = 2 * p.Pout - P[:, p.ny - 1]

    ## ∂P / ∂x = 0 at x = 0, 2
    @. @views P[p.nx, :] = P[p.nx - 1, :]
    @. @views P[1, :] = P[2, :]

    ## ∂P / ∂y = 0 at y = 0 x ∈ [1, 2]
    # @. @views P[nx ÷ 2 + 1: nx, 1] =
    #    P[nx ÷ 2 + 1: nx, 2]
end

function boundary_pressure!(sys::System, p::Parameters)
    boundary_pressure!(sys.P, p)
end

function boundary_velocity!(; f::Function, p::Parameters,
        k::Integer,
        s::AbstractMatrix{T},
        u::AbstractMatrix{T},
        v::AbstractMatrix{T},
        P::AbstractMatrix{T}) where T<:AbstractFloat
    ## u = 0 at x = 0, 2 (walls)
    @. @views u[1, :] = -u[2, :]
    @. @views u[p.nx, :] = -u[p.nx - 1, :]

    ## Darcy at y = 0, 2 (inlet & outlet)
    @. @views v[:, p.ny] = -p.K / p.μ[k] * f(s[:, p.ny]) *
            (P[:, p.ny] - P[:, p.ny - 1]) / p.Δy
    @. @views v[:, 1] =  -p.K / p.μ[k] * f(s[:, 1]) *
            (P[:, 2] - P[:, 1]) / p.Δy
end

function boundary_velocity!(sys::System, p::Parameters)
    ## s is the saturation of gas
    fgas(s) = s^2
    fliq(s) = (1 .- s)^2
    f = [fgas, fliq]

    for k = 1 : 2
        @views boundary_velocity!(f=f[k],
            s=sys.s, u=sys.u[:, :, k], v=sys.v[:, :, k],
            P=sys.P, p=p, k=k)
    end
end
# -------------------------------------------------------- #
