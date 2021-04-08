## Physics:
##
## We simulate 2d flow of ideal gas and C₅H₁₂ through a
## porous medium with the help of Darcy's law.
## 
## v⃗ᵢ = -μᵢ⁻¹ K̂⋅f(sᵢ) ⋅∇P
##
## where K̂, a second-order tensor in the general case, is
## the specific permeability. It depends only on the
## geometry of the mediumm. Here we assume isotropy
## of space, so K is a scalar. μ is the dynamic viscocity.

using DelimitedFiles, PyPlot

## Change index -> [i, j, k]
mutable struct System{T<:AbstractFloat}
    ## index [k, i, j], where k=1,2 is the component: 
    ## 1 - gas,
    ## 2 - liquid;
    ## i, j are x and y spatial coordinates.
    u::Array{T, 3}
    v::Array{T, 3}
    ρ::Array{T, 3}
    ## [i, j] -> x, y
    P::Array{T, 2}
    s::Array{T, 2}
end

Base.copy(sys::System) = System(copy.((sys.u, sys.v, sys.ρ,
                                       sys.P, sys.s))... )

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

function density(x::TaitEoS, P) 
    return x.ρ₀ / (1 - x.C * log10((x.B + P) 
                                   / (x.B + x.P₀)))
end

function density(x::IdealGasEoS, P)
    R = 8.314
    return P * x.M / R / x.T
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

# function densityofLiquid(P::Float64, s::Float64)
#     ρ₀ = 616.18
#     P₀ = 1e5
#     ρ₂ = ρ₀ / (1 - 0.2105 * log10((35e6 + P) / (35e6 + P₀)))
#     ρ̂₂ = (1 - s) * ρ₂
#     return ρ̂₂
# end

## Avoid typos when working with partial derivatives
##
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

## Continuity equation to calculate ρᵢ
## φ ⋅∂ρᵢ / ∂t + div(ρᵢv⃗ᵢ) = 0, where ρᵢ = mᵢ / V
## We will add a right side to this where we inject matter.
function continuity_equation!(; Δt::T, L::T, φ::T,
        ρnext::AbstractMatrix{T},
        ρwork::AbstractMatrix{T},
        uwork::AbstractMatrix{T},
        vwork::AbstractMatrix{T}) where T<:AbstractFloat
    nx, ny = size(ρnext)
    Δx = Δy = L / nx
    for j = 2 : ny - 1
        for i = 2 : nx - 1
            ρc, ρn, ρs, ρw, ρe = cnswe(ρwork, i, j)
            uc, un, us, uw, ue = cnswe(uwork, i, j)
            vc, vn, vs, vw, ve = cnswe(vwork, i, j)
    
            ρnext[i, j] = ρc - Δt / φ *
                (uc * (ρe - ρw) / (2 * Δx) +
                 vc * (ρn - ρs) / (2 * Δy) +
                 ρc * (ue - uw) / (2 * Δx) +
                 ρc * (vn - vs) / (2 * Δy))
        end
    end
end

## Darcy's law
## v⃗ = -μ⁻¹ K̂⋅f(s) ⋅∇P
function darcy!(; f::Function,
        unext::AbstractMatrix{T},
        vnext::AbstractMatrix{T},
        Pnext::AbstractMatrix{T},
            s::AbstractMatrix{T},
        K::T, μ::T, L::T) where T<:AbstractFloat
    nx, ny = size(unext)
    Δx = Δy = L / nx
    for j = 2 : ny - 1
        for i = 2 : nx - 1
            ## можно через (ρwork + ρnext) / 2
            Pc, Pn, Ps, Pw, Pe = cnswe(Pnext, i, j)

            unext[i, j] = -K / μ * f(s[i, j]) *
                (Pe - Pw) / (2 * Δx)
            vnext[i, j] = -K / μ * f(s[i, j]) *
                (Pn - Ps) / (2 * Δy)
        end
    end
end

function boundaryDensity!(ρ::AbstractArray{T, 3},
        P::AbstractMatrix{T}) where T<:AbstractFloat
    nx, ny = size(ρ)

    @. @views ρ[1, :, 1]  = density(ideal_gas, P[:, 1])
    @. @views ρ[1, :, ny] = density(ideal_gas, P[:, ny])
    @. @views ρ[1, 1, :]  = density(ideal_gas, P[1, :])
    @. @views ρ[1, nx, :] = density(ideal_gas, P[nx, :])

    @. @views ρ[2, :, 1]  = density(tait_C₅H₁₂, P[:, 1])
    @. @views ρ[2, :, ny] = density(tait_C₅H₁₂, P[:, ny])
    @. @views ρ[2, 1, :]  = density(tait_C₅H₁₂, P[1, :])
    @. @views ρ[2, nx, :] = density(tait_C₅H₁₂, P[nx, :])
end

function boundarySaturation!(ψ::T, M₁::T, M₂::T,
        ρ::AbstractArray{T, 3},
        s::AbstractMatrix{T}) where T<:AbstractFloat
    nx, ny = size(ρ)
    ## Inlet
    ## We derive the saturation from the molar composition and
    ## the densities.
    @. @views s[1: nx ÷ 2, 1] = (ρ[1, 1: nx ÷ 2, 1] /
            ρ[2, 1: nx ÷ 2, 1] * M₁ / M₂ / ψ + 1)^(-1)
end

function initialSaturation!(ψ₀::T, M₁::T, M₂::T,
        ρ::AbstractArray{T, 3},
        s::AbstractMatrix{T}) where T<:AbstractFloat
    @. @views s[:, :] = (ρ[1, :, :] /
            ρ[2, :, :] * M₁ / M₂ / ψ₀ + 1)^(-1)
end

function boundaryPressure!(P::Array{T, 2},
        Pin::T, Pout::T) where T<:AbstractFloat
    nx, ny = size(P)
    ## P = Pin at y = 0 and x ∈ [0, 1]
    @. @views P[:, 1] = 2 * Pin - P[:, 2]

    ## P = Pout at y = 2
    @. @views P[:, ny] = 2 * Pout - P[:, ny - 1]

    ## ∂P / ∂x = 0 at x = 0, 2
    @. @views P[nx, :] = P[nx - 1, :]
    @. @views P[1, :] = P[2, :]

    ## ∂P / ∂y = 0 at y = 0 x ∈ [1, 2]
    # @. @views P[nx ÷ 2 + 1: nx, 1] =
    #    P[nx ÷ 2 + 1: nx, 2]
end

function boundaryVelocity!(; K::T, μ::T, f::Function, L::T,
        s::AbstractMatrix{T},
        u::AbstractMatrix{T},
        v::AbstractMatrix{T},
        P::AbstractMatrix{T}) where T<:AbstractFloat
    nx, ny = size(u)
    Δy = L / nx
    ## u = 0 at x = 0, 2 (walls)
    @. @views u[1, :] = -u[2, :]
    @. @views u[nx, :] = -u[nx - 1, :]

    ## v = 0 at y = 0 and x ∈ [1, 2]
    # @. @views v[1, nx ÷ 2 + 1: nx] =
    #     -v[2, nx ÷ 2 + 1: nx]

    ## Darcy 
    @. @views v[:, ny] = -K / μ * f(s[:, ny]) *
        (P[:, ny] - P[:, ny - 1]) / Δy
    @. @views v[:, 1] =  -K / μ * f(s[:, 1]) *
        (P[:, 2] - P[:, 1]) / Δy
end

function update!(; ρnext, unext, vnext, Pnext, snext,
        ρwork, uwork, vwork,
        Δt, L, Pin, Pout, φ, K, μ, M, ψ)
    nx, ny = size(ρnext)
    Δx = Δy = L / nx

    for k = 1 : 2
        @views continuity_equation!(Δt=Δt, L=L, φ=φ,
            ρnext=ρnext[k, :, :], ρwork=ρwork[k, :, :],
            uwork=uwork[k, :, :], vwork=vwork[k, :, :])
    end

    for index in CartesianIndices(Pnext)
        i, j = Tuple(index)
        Pnext[i, j], snext[i, j] =
            findPressure(ρ̂₁=ρnext[1, i, j],
                         ρ̂₂=ρnext[2, i, j], V=Δx * Δy)
    end

    boundaryPressure!(Pnext, Pin, Pout)
    boundaryDensity!(ρnext, Pnext) 
    boundarySaturation!(ψ, M[1], M[2], ρnext, snext) 
    
    ## s is the saturation of gas
    fgas(s) = s^2
    fliq(s) = (1 .- s)^2
    f = [fgas, fliq]

    for k = 1 : 2
        @views darcy!(f=f[k], unext=unext[k, :, :],
               vnext=vnext[k, :, :], Pnext=Pnext, s=snext,
               K = K, μ=μ[k], L=L)

        @views boundaryVelocity!(f=f[k], K=K, μ=μ[k], L=L,
            s=snext, u=unext[k, :, :], v=vnext[k, :, :],
            P=Pnext)
    end
end

function filtration!(; L, nx, ny, Δt, nsteps, K, φ, μ,
        Pin, Pout, P₀, M, ψ, ψ₀)
    Δx = Δy = L / nx

    P = fill(P₀, nx, ny)
    boundaryPressure!(P, Pin, Pout)

    s = zeros(nx, ny)
    ρ = zeros(2, nx, ny)

    @views ρ[2, :, :] .= density(tait_C₅H₁₂, P₀)
    @views ρ[1, :, :] .= density(ideal_gas, P₀)

    initialSaturation!(ψ₀, M[1], M[2], ρ, s)

    u = zeros(2, nx, ny)
    v = zeros(2, nx, ny)
    
    ## s is the saturation of gas
    fgas(s) = s^2
    fliq(s) = (1 .- s)^2
    f = [fgas, fliq]

    for k = 1 : 2
        @views boundaryVelocity!(f=f[k], K=K, μ=μ[k], P=P,
                    L=L, u=u[k, :, :], v=v[k, :, :], s=s)
    end

    syswork = System(u, v, ρ, P, s)
    sysnext = copy(syswork)
    
    for t = 1 : nsteps
        Parameters = (
                      ρnext = sysnext.ρ,
                      unext = sysnext.u,
                      vnext = sysnext.v,
                      Pnext = sysnext.P,
                      snext = sysnext.s,
                      ρwork = syswork.ρ,
                      uwork = syswork.u,
                      vwork = syswork.v,
                      Δt = Δt,
                      L = L,
                      Pin = Pin,
                      Pout = Pout,
                      φ = φ,
                      K = K,
                      μ = μ,
                      M = M,
                      ψ = ψ
                     )
        update!(; Parameters...)

        ## Swap next and work and continue iteration
        sysnext, syswork = syswork, sysnext
    end
# 
#     ## If in the end syswork is not the original sys,
#     ## copy the most recent data into sys
#     if syswork !== sys
#         sys = syswork
#     end
    return syswork
end

function plot(u::Array{T, 3}, v::Array{T, 3},
        P::Array{T, 2}) where T<:AbstractFloat
    nx, ny = size(u[1, :, :])

    fig = PyPlot.figure(figsize=(10, 10))
    ax = PyPlot.axes()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    rgx = range(0, 2, length=nx)
    rgy = range(0, 2, length=ny)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))
    
    ## Velocity vector field
    ## scale = 100.0
    quiver(x, y, u[1, :, :]', v[1, :, :]', color="r")
    quiver(x, y, u[2, :, :]', v[2, :, :]', color="b")
    
    ## Pressure contour map
    pos = ax.contourf(X, Y, P', alpha=0.5,
        cmap=matplotlib.cm.viridis)
    fig.colorbar(pos, ax=ax)
    cp = contour(X, Y, P', cmap=matplotlib.cm.viridis)

    PyPlot.title("2 Phase Filtration: Liquid and Ideal Gas")
    
    savefig("img/2phase-filtration.png", dpi=200)
    savefig("img/2phase-filtration.svg")
end

let
    duration = 2000
    Δt = 0.000001
    nsteps = round(Int64, duration / Δt)
    nsteps = 10
    
    ## Space grid [0, 2] × [0, 2]
    ## Change later to 2 (how ?)
    L = 2.0
    nx = ny = 40

    ## Porosity
    φ = 0.7
    
    ## Specific permeability
    K = 1e-12
    
    ## Dynamic viscocity
    μ = [1.8e-5, 2.14e-4]
    
    ## Pressures at top and bottom, and initial pressure
    Pin = 1e6
    Pout = 1e5
    P₀ = Pout

    ## Molar masses
    M = [0.028, 0.07215]

    ## Molar Composition on the Inlet
    ## ψ = ν₁ / ν₂ = (m₁ / M₁) / (m₂ / M₂)
    ψ = 0.3

    ## Initial Molar Composition
    ψ₀ = ψ

    Parameters = (
                  Δt = Δt,
                  nsteps = nsteps,
                  nx = nx,
                  ny = ny,
                  L = L,
                  φ = φ,
                  K = K,
                  μ = μ,
                  Pin = Pin,
                  Pout = Pout,
                  P₀ = P₀,
                  M = M,
                  ψ = ψ,
                  ψ₀ = ψ₀
                 )

    sys = filtration!(; Parameters...)
    plot(sys.u, sys.v, sys.P)
    # writedlm("pressure.txt", sys.P, ' ')
    # writedlm("u1.txt", sys.u[1, :, :], ' ')
    # writedlm("u2.txt", sys.u[2, :, :], ' ')
    # writedlm("density1.txt", sys.ρ[1, :, :], ' ')
    # writedlm("density2.txt", sys.ρ[2, :, :], ' ')
end
