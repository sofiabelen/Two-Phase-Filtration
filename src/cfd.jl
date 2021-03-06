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

# const fgas(s::AbstractFloat)::AbstractFloat = s^2
# const fliq(s::AbstractFloat)::AbstractFloat = (1 .- s)^2
# const fresisting = (fgas, fliq)::Tuple{Function, Vararg{Function}}

function f_phase(s::AbstractFloat)
    return s^2
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
    @debug ρnext ρwork
end

function continuity_equation!(syswork::System, 
        sysnext::System, p::Parameters)
    for k = 1 : 2
        @debug "Component $k continuity eqn"
        @views continuity_equation!(ρnext=sysnext.ρ[:, :, k], 
                                    ρwork=syswork.ρ[:, :, k],
                                    uwork=syswork.u[:, :, k],
                                    vwork=syswork.v[:, :, k],
                                    p=p)
    end
    @debug "Saturation " sysnext.s
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
    ## Same problem with type-stability as boundary_velocity
    ## https://discourse.julialang.org/t/how-to-make-type-stable-passing-an-element-of-a-vector-of-functions/60068
    # fgas(s) = s^2
    # fliq(s) = (1 .- s)^2
    # f = (fgas, fliq)

    f(s::AbstractFloat)::AbstractFloat = s^2

    for k = 1 : 2
      @views darcy!(f=f,
                    u=sys.u[:, :, k],
                    v=sys.v[:, :, k], P=sys.P,
                    s=sys.s[:, :, k], p=p, k=k)
    end

    ## This is type stable
    ## (but I'm still getting runtime dispatch in profiler)
    # fgas(s) = s^2
    # fliq(s) = (1 .- s)^2

    # @views darcy!(f=fgas, u=sys.u[:, :, 1],
    #               v=sys.v[:, :, 1], P=sys.P,
    #               s=sys.s, p=p, k=1)

    # @views darcy!(f=fliq, u=sys.u[:, :, 2],
    #               v=sys.v[:, :, 2], P=sys.P,
    #               s=sys.s, p=p, k=2)
end
# -------------------------------------------------------- #



# ---------------------- BC & IC ------------------------- #

function boundary_density!(ρ::AbstractArray{T, 3},
        P::AbstractMatrix{T}, s::AbstractArray{T, 3},
        p::Parameters) where T<:AbstractFloat
    @. @views ρ[:, 1, 1]  = density(ideal_gas, P[:, 1]) *
                                                 s[:, 1, 1]
    @. @views ρ[:, p.ny, 1] = density(ideal_gas, P[:, p.ny])*
                                                 s[:, p.ny, 1]
    @. @views ρ[1, :, 1]  = density(ideal_gas, P[1, :]) *
                                                 s[1, :, 1]
    @. @views ρ[p.nx, :, 1] = density(ideal_gas, P[p.nx, :])*
                                                 s[p.nx, :, 1]  

    @. @views ρ[:, 1, 2]  = density(tait_C₅H₁₂, P[:, 1]) *
                                                 s[p.nx, :, 2]
    @. @views ρ[:, p.ny, 2] = density(tait_C₅H₁₂, P[:,p.ny])*
                                                 s[p.nx, :, 2]
    @. @views ρ[1, :, 2]  = density(tait_C₅H₁₂, P[1, :]) *
                                                 s[p.nx, :, 2]
    @. @views ρ[p.nx, :, 2] = density(tait_C₅H₁₂, P[p.nx,:])*
                                                 s[p.nx, :, 2]
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
        s::AbstractArray{T, 3},
        p::Parameters) where T<:AbstractFloat

    ## This is type unstable.
    @. @views s[1: p.nx ÷ 2, 1, 1] = (ρ̂₁[1: p.nx ÷ 2, 1] /
                                   ρ̂₂[1: p.nx ÷ 2, 1] *
                         p.M[2] * (1 - p.ψ) / p.M[1] /
                         p.ψ + 1)^(-1)
    @. @views s[:, ny, 1] = s[:, ny - 1, 1]

    @. @views s[1: p.nx ÷ 2, 1, 2] = 1 - s[1: p.nx ÷ 2, 1, 1]
    @. @views s[:, ny, 1] = 1 - s[:, ny, 1]
end

function boundary_saturation!(sys::System, p::Parameters)
    ρ̂₁ = density.(ideal_gas, sys.P)
    ρ̂₂ = density.(tait_C₅H₁₂, sys.P)
    boundary_saturation!(ρ̂₁, ρ̂₂, sys.s, p)
end

function initial_saturation!(
        ρ̂₁::AbstractMatrix{T},
        ρ̂₂::AbstractMatrix{T},
        s::AbstractArray{T, 3}, p) where T<:AbstractFloat

    @. @views s[:, :, 1] = (ρ̂₁[:, :] /
                         ρ̂₂[:, :] * p.M[2] * (1 - p.ψ₀) /
                         p.M[1] / p.ψ₀ + 1)^(-1)

    @. @views s[:, :, 2] = 1 - s[:, :, 1]
end

function initial_saturation!(sys::System, p::Parameters)
    ρ̂₁ = density.(ideal_gas, sys.P)
    ρ̂₂ = density.(tait_C₅H₁₂, sys.P)
    initial_saturation!(ρ̂₁, ρ̂₂, sys.s, p)
end

function boundary_pressure!(P::Array{T, 2},
        p::Parameters) where T<:AbstractFloat
    ## P = Pin at y = 0 && x ∈ [0, 1]
    @. @views P[1 : nx ÷ 2, 1] = 2 * p.Pin - P[1 : nx ÷ 2, 2]

    ## P = Pout at y = 2
    @. @views P[:, p.ny] = 2 * p.Pout - P[:, p.ny - 1]

    ## ∂P / ∂x = 0 at x = 0, 2
    @. @views P[p.nx, :] = P[p.nx - 1, :]
    @. @views P[1, :] = P[2, :]

    ## ∂P / ∂y = 0 at y = 0  && x ∈ [1, 2]
    @. @views P[nx ÷ 2 + 1: nx, 1] =
       P[nx ÷ 2 + 1: nx, 2]
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
    ## v = 0 at y = 0 && x ∈ [1, 2]
    @. @views v[p.nx ÷ 2 + 1: p.nx, 1] = 0

    ## Darcy at y = 0 && x ∈ [0, 1], y = 2 (inlet & outlet)
    @. @views v[1: nx ÷ 2, 1] =  -p.K / p.μ[k] *
     f(s[1: nx ÷ 2, 1]) * (-3 * P[1 : nx ÷ 2, 1] +
     4 * P[1 : nx ÷ 2, 2] - P[1 : nx ÷ 2, 3]) / (2 * p.Δy)

    @. @views v[:, p.ny] = -p.K / p.μ[k] * f(s[:, p.ny]) *
    (-3 * P[:, ny] + 4 * P[:, ny - 1] - P[:, ny - 2]) /
    (-2 * p.Δy)
end

function boundary_velocity!(sys::System, p::Parameters)
    ## Type unstable because of 'f=f[k]'
    ## This answer doesn't work because of views
    ## https://discourse.julialang.org/t/how-to-make-type-stable-passing-an-element-of-a-vector-of-functions/60068

    f(s::AbstractFloat)::AbstractFloat = s^2

    for k = 1 : 2
        @views boundary_velocity!(f=f, s=sys.s[:, :, k],
            u=sys.u[:, :, k], v=sys.v[:, :, k],
            P=sys.P, p=p, k=k)
    end

    # fgas(s) = s^2
    # fliq(s) = (1 .- s)^2

    ## This is type stable
    # @views boundary_velocity!(f=fgas,
    #     s=sys.s, u=sys.u[:, :, 1], v=sys.v[:, :, 1],
    #     P=sys.P, p=p, k=1)

    # @views boundary_velocity!(f=fliq,
    #     s=sys.s, u=sys.u[:, :, 2], v=sys.v[:, :, 2],
    #     P=sys.P, p=p, k=2)
end
# -------------------------------------------------------- #
