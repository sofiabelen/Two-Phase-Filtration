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

function f_phase(s::AbstractFloat)::AbstractFloat
    return s^2
end


# ------------------ Continuity Equation ----------------- #

## Continuity equation to calculate ρᵢ
## φ ⋅∂ρᵢ / ∂t + div(ρᵢv⃗ᵢ) = 0, where ρᵢ = mᵢ / V
## (We can add a right side to this if we inject matter.)
function continuity_equation!(;
        ρnext::AbstractMatrix{T},
        ρwork::AbstractMatrix{T},
        F::AbstractMatrix{T},
        p::Parameters) where T<:AbstractFloat
    for j = 2 : p.ny - 1
        for i = 2 : p.nx - 1
            ρc, ρn, ρs, ρw, ρe = cnswe(ρwork, i, j)
            ρnext[i, j] = ρc + p.Δt  * F[i, j]
        end
    end
    @debug ρnext ρwork
end

function continuity_equation_fill_F!(;
        ρwork::AbstractMatrix{T},
        uwork::AbstractMatrix{T},
        vwork::AbstractMatrix{T},
        F::AbstractMatrix{T},
        p::Parameters) where T<:AbstractFloat
    for j = 2 : p.ny - 1
        for i = 2 : p.nx - 1
            ρc, ρn, ρs, ρw, ρe = cnswe(ρwork, i, j)
            uc, un, us, uw, ue = cnswe(uwork, i, j)
            vc, vn, vs, vw, ve = cnswe(vwork, i, j)
    
            F[i, j] = - 1 / p.φ *
                (uc * (ρe - ρw) / (2 * p.Δx) +
                 vc * (ρn - ρs) / (2 * p.Δy) +
                 ρc * (ue - uw) / (2 * p.Δx) +
                 ρc * (vn - vs) / (2 * p.Δy))
        end
    end
    @debug ρnext ρwork
end

## Predictor
function continuity_equation!(syswork::System, 
        sysnext::System, p::Parameters)
    for k = 1 : 2
        @debug "Component $k continuity eqn"
        @views continuity_equation_fill_F!(
                                    ρwork=syswork.ρ[:, :, k],
                                    uwork=syswork.u[:, :, k],
                                    vwork=syswork.v[:, :, k],
                                    F=syswork.F[:, :, k], p=p)
        @views continuity_equation!(ρnext=sysnext.ρ[:, :, k], 
                                    ρwork=syswork.ρ[:, :, k],
                                    F=syswork.F[:, :, k], p=p)
    end
    @debug "Saturation " sysnext.s
end

## Corrector
function continuity_equation!(syswork::System,
        sysnext_pred::System, sysnext::System,
        p::Parameters)
    for k = 1 : 2
        @debug "Component $k continuity eqn"
        @views continuity_equation_fill_F!(
                                    ρwork=sysnext_pred.ρ[:, :, k],
                                    uwork=sysnext_pred.u[:, :, k],
                                    vwork=sysnext_pred.v[:, :, k],
                                    F=sysnext_pred.F[:, :, k], p=p)
        @views continuity_equation!(ρnext=sysnext.ρ[:, :, k], 
                                    ρwork=syswork.ρ[:, :, k],
                                    F=(syswork.F[:, :, k] .+
                                       sysnext_pred.F[:, :, k]) ./2,
                                    p=p)
    end
    @debug "Saturation " sysnext.s
end
# -------------------------------------------------------- #


# ---------------------- Darcy --------------------------- #

## Darcy's law
## v⃗ = -μ⁻¹ K̂⋅f(s) ⋅∇P

function darcy!(;
        u::AbstractMatrix{T}, v::AbstractMatrix{T},
        P::AbstractMatrix{T}, s::AbstractMatrix{T},
        p::Parameters, k::Integer) where T<:AbstractFloat
    for j = 2 : p.ny - 1
        for i = 2 : p.nx - 1
            ## можно через (ρwork + ρ) / 2 -> 2-го порядка
            Pc, Pn, Ps, Pw, Pe = cnswe(P, i, j)

            u[i, j] = -p.K / p.μ[k] * f_phase(s[i, j]) *
                (Pe - Pw) / (2 * p.Δx)
                v[i, j] = -p.K / p.μ[k] * f_phase(s[i, j]) *
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

    # f(s::AbstractFloat)::AbstractFloat = s^2

    for k = 1 : 2
      @views darcy!(u=sys.u[:, :, k],
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

function boundary_saturation!(
        ρ̂₁::AbstractMatrix{T},
        ρ̂₂::AbstractMatrix{T},
        s::AbstractArray{T, 3},
        p::Parameters) where T<:AbstractFloat
## Inlet
## We derive the saturation from the molar composition
## ψ = ν₁ / ν₂ and the densities ρ̂₁ and ρ̂₂.
## s_boundary = (s_1 + s_2) / 2
    ρ̂₁_bound::Array{AbstractFloat} = ρ̂₁[:, 1]
    ρ̂₂_bound::Array{AbstractFloat} = ρ̂₁[:, 1]
    s_bound::Array{AbstractFloat}  = ρ̂₁[:, 1]

    @. @views ρ̂₁_bound = (ρ̂₁[:, 1]
                + ρ̂₁[:, 2]) / 2
    @. @views ρ̂₂_bound = (ρ̂₂[:, 1]
                + ρ̂₂[:, 2]) / 2
    @. @views s_bound  = (ρ̂₁_bound / ρ̂₂_bound *
               p.M[2] * (1 - p.ψ) / p.M[1] / p.ψ + 1)^(-1)
    @. @views s[:, 1, 1] = 2 * s_bound -
               s[:, 1, 1]

## Outlet
## ∂s / ∂n⃗ = 0
    @. @views s[:, ny, 1] = s[:, ny - 1, 1]

    @. @views s[:, 1, 2] = 1 - s[:, 1, 1]
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
    @. @views P[:, 1] = 2 * p.Pin - P[:, 2]

    ## P = Pout at y = 2
    @. @views P[:, p.ny] = 2 * p.Pout - P[:, p.ny - 1]

    ## ∂P / ∂x = 0 at x = 0, 2
    @. @views P[p.nx, :] = P[p.nx - 1, :]
    @. @views P[1, :] = P[2, :]

#     ## ∂P / ∂y = 0 at y = 0  && x ∈ [1, 2]
#     @. @views P[nx ÷ 2 + 1: nx, 1] =
#        P[nx ÷ 2 + 1: nx, 2]
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
#     ## v = 0 at y = 0 && x ∈ [1, 2]
#     @. @views v[p.nx ÷ 2 + 1: p.nx, 1] = 0

    ## Darcy at y = 0 && x ∈ [0, 1], y = 2 (inlet & outlet)
    @. @views v[:, 1] =  -p.K / p.μ[k] *
     f_phase(s[:, 1]) * (-3 * P[:, 1] +
     4 * P[:, 2] - P[:, 3]) / (2 * p.Δy)

    @. @views v[:, p.ny] = -p.K / p.μ[k] * f_phase(s[:, p.ny]) *
    (-3 * P[:, ny] + 4 * P[:, ny - 1] - P[:, ny - 2]) /
    (-2 * p.Δy)
end

function boundary_velocity!(sys::System, p::Parameters)
    ## Type unstable because of 'f=f[k]'
    ## This answer doesn't work because of views
    ## https://discourse.julialang.org/t/how-to-make-type-stable-passing-an-element-of-a-vector-of-functions/60068

    # f(s::AbstractFloat)::AbstractFloat = s^2

    for k = 1 : 2
        @views boundary_velocity!(s=sys.s[:, :, k],
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

# --------------------- New BC --------------------------- #
function apply_bc!(; grid_ghost::AbstractMatrix{T},
        grid_inside::AbstractMatrix{T}, a::AbstractMatrix{T},
        b::AbstractMatrix{T}, rhs::AbstractMatrix{T},
        p::Parameters, Δh::T) where T<:AbstractFloat
    @views @. grid_ghost = (2 * Δh * rhs + grid_inside *
                                   (2 * b - a * Δh)) /
                                  (a * Δh + 2 * b)
    # println("Denominator: ", (a * p.Δx + 2 * b))
    # println("ghost ", grid_ghost)
    # println("")
end

function apply_bc!(grid::Array{T, 3}, bc::BoundaryCondition{T},
        rhs::Array{T, 3}, p::Parameters) where T<:AbstractFloat
    ## left
    # println("left")
    @views apply_bc!(grid_ghost  = grid[1, :, :],
              grid_inside = grid[2, :, :],
              a   = bc.a[:, 1, :],
              b   = bc.b[:, 1, :],
              rhs =  rhs[:, 1, :], p=p, Δh = p.Δx)

    ## right
    # println("right")
    @views apply_bc!(grid_ghost  = grid[nx,     :, :],
              grid_inside = grid[nx - 1, :, :],
              a   = bc.a[:, 2, :],
              b   = bc.b[:, 2, :],
              rhs =  rhs[:, 2, :], p=p, Δh = p.Δx)

    ## top
    # println("top")
    @views apply_bc!(grid_ghost  = grid[:, ny,     :],
              grid_inside = grid[:, ny - 1, :],
              a   = bc.a[:, 3, :],
              b   = bc.b[:, 3, :],
              rhs =  rhs[:, 3, :], p=p, Δh = p.Δy)

    ## bottom
    # println("bottom")
    @views apply_bc!(grid_ghost  = grid[:, 1, :],
              grid_inside = grid[:, 2, :],
              a   = bc.a[:, 4, :],
              b   = bc.b[:, 4, :],
              rhs =  rhs[:, 4, :], p=p, Δh = p.Δy)
end

function update_RHS_saturation!(; s_RHS::AbstractMatrix{T},
        P::Array{T, 1}, p::Parameters) where T<:AbstractFloat
    @views ρ̂₁_bound    = density.(ideal_gas,  P)
    @views ρ̂₂_bound    = density.(tait_C₅H₁₂, P)
    @views s_RHS[:, 1] = saturation.(ρ̂₁_bound, ρ̂₂_bound, p.ψ,
                                         p.M[1], p.M[2])
    @views @. s_RHS[:, 2] = 1 - s_RHS[:, 1]
end

function update_RHS_saturation!(s_RHS::Array{T, 3}, P::Array{T, 2},
        p::Parameters) where T<:AbstractFloat
    ## bottom inlet
    @views update_RHS_saturation!(s_RHS = s_RHS[p.inlet, 4, :],
                           P = fill(p.Pin, last(collect(p.inlet))),
                           p = p)
    ## everywhere else we don't have to do anything because
    ## ∂s / ∂n = 0
end

function get_bound(h::AbstractArray{T, 2},
        p::Parameters) where T<:AbstractFloat
    @views h_left  = (h[1, :]    + h[2, :])        / 2
    @views h_right = (h[p.nx, :] + h[p.nx - 1, :]) / 2
    @views h_top   = (h[:, p.ny] + h[:, p.ny - 1]) / 2
    @views h_bott  = (h[:, 1]    + h[:, 2])        / 2
    return hcat(h_left, h_right, h_top, h_bott)
end

function update_RHS_density!(ρ_RHS::Array{T, 3},
        s::Array{T, 3}, P::Array{T, 2},
        p::Parameters) where T<:AbstractFloat
    @views P_bound  = get_bound(P, p)
    @views s₁_bound = get_bound(s[:, :, 1], p)
    @views s₂_bound = 1 .- s₁_bound

    @views ρ̂₁_bound = density.(ideal_gas,  P_bound)
    @views ρ̂₂_bound = density.(tait_C₅H₁₂, P_bound)
    @views @. ρ_RHS[:, :, 1] = s₁_bound * ρ̂₁_bound
    @views @. ρ_RHS[:, :, 2] = s₂_bound * ρ̂₂_bound
end

## We call this only once during init
function update_RHS_pressure!(P_RHS::Array{T, 3},
        p::Parameters) where T<:AbstractFloat
    ## bottom inlet
    @views @. P_RHS[p.inlet, 4, 1] = p.Pin

    ## top outlet
    @views @. P_RHS[:, 3, 1]       = p.Pout

    ## at walls ∂P / ∂n = 0 so we don't have to change anything,
    ## since P_RHS = 0
end

function update_RHS_velocity!(v_RHS::Array{T, 3}, 
        P::Array{T, 2}, s::Array{T, 3},
        p::Parameters) where T<:AbstractFloat
    @views s₁_bound = get_bound(s[:, :, 1], p)
    @views s₂_bound = 1 .- s₁_bound
    @views s_bound = reshape(hcat(s₁_bound, s₂_bound),
                      (max(p.nx, p.ny), 4, 2))

    ## bottom inlet
    for k = 1 : 2
        @. @views v_RHS[p.inlet, 4, k] =  -p.K / p.μ[k] *
                                  f_phase(s_bound[p.inlet, 4, k]) *
                                  (P[p.inlet, 2] -
                                   P[p.inlet, 1]) / p.Δy

        ## top outlet
        @. @views v_RHS[:, 3, k]       =  -p.K / p.μ[k] *
                                    f_phase(s_bound[:, 3, k]) *
                                    (P[:, ny] -
                                     P[:, ny - 1]) / p.Δy
    end
end

function apply_bc_pres_sat_den!(sys::System, p::Parameters)
    ## ∂P / ∂n = 0 at walls
    ## P = Pin     at inlet
    ## P = Pout    at outlet
    P_tmp = zeros(p.nx, p.ny, 2)
    P_tmp[:, :, 1] = sys.P
    apply_bc!(P_tmp, p.bc.P, sys.bc_RHS.P, p)
    sys.P = P_tmp[:, :, 1]

    ## s_inlet = s_inlet(ρ̂₁(T, P), ρ̂₂(P, T), P, T)
    update_RHS_saturation!(sys.bc_RHS.s, sys.P, p)
    apply_bc!(sys.s, p.bc.s, sys.bc_RHS.s, p)

    ## ρ_k_bound = ρ̂_k_bound(T, P) * s_k_bound(T, P)
    update_RHS_density!(sys.bc_RHS.ρ, sys.s, sys.P, p)
    apply_bc!(sys.ρ, p.bc.ρ, sys.bc_RHS.ρ, p)
end

function apply_bc_velocity!(sys::System, p::Parameters)
    ## Darcy at inlet and outlet
    update_RHS_velocity!(sys.bc_RHS.v, sys.P, sys.s, p)
    apply_bc!(sys.u, p.bc.u, sys.bc_RHS.u, p)
    apply_bc!(sys.v, p.bc.v, sys.bc_RHS.v, p)
end

function initial_conditions!(sys::System, p::Parameters)
    ## ∂P / ∂n = 0 at walls
    ## P = Pin     at inlet
    ## P = Pout    at outlet
    update_RHS_pressure!(sys.bc_RHS.P, p)
    P_tmp = zeros(p.nx, p.ny, 2)
    P_tmp[:, :, 1] = sys.P
    apply_bc!(P_tmp, p.bc.P, sys.bc_RHS.P, p)
    sys.P = P_tmp[:, :, 1]

    initial_saturation!(sys, p)
    # println("Initial saturation")
    # println(sys.s)

    ## s_inlet = s_inlet(ρ̂₁(T, P), ρ̂₂(P, T), P, T)
    update_RHS_saturation!(sys.bc_RHS.s, sys.P, p)
    apply_bc!(sys.s, p.bc.s, sys.bc_RHS.s, p)

    # println("Saturation after bc:")
    # println(sys.s)

    # println("Saturation RHS bottom")
    # println(sys.bc_RHS.s[:, 4, 1])

    # println("Saturatin inlet ghost bottom")
    # println(sys.s[p.inlet, 1, 1])

    # println("Saturatin inlet inside bottom")
    # println(sys.s[p.inlet, 2, 1])

    density!(sys)
    update_RHS_density!(sys.bc_RHS.ρ, sys.s, sys.P, p)
    apply_bc!(sys.ρ, p.bc.ρ, sys.bc_RHS.ρ, p)

    ## Darcy at inlet and outlet
    update_RHS_velocity!(sys.bc_RHS.v, sys.P, sys.s, p)
    # println("RHS_velocity", sys.bc_RHS.v)
    apply_bc!(sys.u, p.bc.u, sys.bc_RHS.u, p)
    apply_bc!(sys.v, p.bc.v, sys.bc_RHS.v, p)
end

function initial_RHS(p::Parameters)
    u_RHS = zeros(max(p.nx, p.ny), 4, 2)
    v_RHS = zeros(max(p.nx, p.ny), 4, 2)
    ρ_RHS = zeros(max(p.nx, p.ny), 4, 2)
    s_RHS = zeros(max(p.nx, p.ny), 4, 2)
    P_RHS = zeros(max(p.nx, p.ny), 4, 2)
    return BoundaryRHS(u_RHS, v_RHS, ρ_RHS, s_RHS, P_RHS)
end
# -------------------------------------------------------- #
