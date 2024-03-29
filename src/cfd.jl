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

## Fluxes -> staggered grid
function nswe(Φu::AbstractMatrix{T},
        Φv::AbstractMatrix{T},
        i::Int, j::Int) where T<:AbstractFloat
    s = Φv[i, j - 1]
    w = Φu[i - 1, j]
    n = Φv[i, j]
    e = Φu[i, j]
    return n, s, w, e
end

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
            ρnext[i, j] = ρwork[i, j] + p.Δt * F[i, j]
        end
    end
    @debug ρnext ρwork
end

function continuity_equation_fill_F!(;
        Φu::AbstractMatrix{T},
        Φv::AbstractMatrix{T},
        F::AbstractMatrix{T},
        p::Parameters) where T<:AbstractFloat
    for j = 2 : p.ny - 1
        for i = 2 : p.nx - 1
            Φn, Φs, Φw, Φe = nswe(Φu, Φv, i, j)

            F[i, j] = - 1 / p.φ * ((Φe - Φw) / p.Δx +
                                   (Φn - Φs) / p.Δy)
        end
    end
end

# -------------------- Predictor ------------------------- #
function continuity_equation!(syswork::System, 
        sysnext::System, p::Parameters)
    for k = 1 : 2
        @debug "Component $k continuity eqn"
        @views continuity_equation_fill_F!(
                                    Φu=syswork.Φu[:, :, k],
                                    Φv=syswork.Φv[:, :, k],
                                    F=syswork.F[:, :, k], p=p)

        @views continuity_equation!(ρnext=sysnext.ρ[:, :, k], 
                                    ρwork=syswork.ρ[:, :, k],
                                    F=syswork.F[:, :, k], p=p)
    end
    @debug "Saturation " sysnext.s
end
# -------------------------------------------------------- #


# -------------------- Corrector ------------------------- #
function continuity_equation!(syswork::System,
        sysnext_pred::System, sysnext::System,
        p::Parameters)
    for k = 1 : 2
        @views continuity_equation_fill_F!(
                                    Φu=sysnext_pred.Φu[:, :, k],
                                    Φv=sysnext_pred.Φv[:, :, k],
                                    F=sysnext_pred.F[:, :, k], p=p)

        ## (F + F̃) / 2
        @views F_avg = (syswork.F[:, :, k] .+
                        sysnext_pred.F[:, :, k]) ./2

        @views continuity_equation!(ρnext=sysnext.ρ[:, :, k], 
                                    ρwork=syswork.ρ[:, :, k],
                                    F=F_avg, p=p)
    end
end
# -------------------------------------------------------- #


# ---------------------- Darcy --------------------------- #

## Darcy's law
## v⃗ = -μ⁻¹ K̂⋅f(s) ⋅∇P
function darcy!(;
        u::AbstractMatrix{T}, v::AbstractMatrix{T},
        P::AbstractMatrix{T}, s::AbstractMatrix{T},
        p::Parameters, k::Integer) where T<:AbstractFloat
    ## We iterate through the whole u and v arrays, so it turns
    ## out, that we don't need BC on velocities thanks to the
    ## staggered grid :)
    for j = 1 : p.ny
        for i = 1 : p.nx - 1
            Pw = P[i, j]
            Pe = P[i + 1, j]
            u[i, j] = -p.K / p.μ[k] * f_phase(s[i, j]) *
                (Pe - Pw) / p.Δx

        end
    end

    for j = 1 : p.ny - 1
        for i = 1 : p.nx
            Ps = P[i, j]
            Pn = P[i, j + 1]
            v[i, j] = -p.K / p.μ[k] * f_phase(s[i, j]) *
                (Pn - Ps) / p.Δy
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

function apply_bc!(; grid_ghost::AbstractMatrix{T},
        grid_inside::AbstractMatrix{T}, a::AbstractMatrix{T},
        b::AbstractMatrix{T}, rhs::AbstractMatrix{T},
        p::Parameters, Δh::T) where T<:AbstractFloat
    @views @. grid_ghost = (2 * Δh * rhs + grid_inside *
                                   (2 * b - a * Δh)) /
                                  (a * Δh + 2 * b)
end

function apply_bc!(grid::Array{T, 3}, bc::BoundaryCondition{T},
        rhs::Array{T, 3}, p::Parameters) where T<:AbstractFloat
    ## left
    @views apply_bc!(grid_ghost  = grid[1, :, :],
              grid_inside = grid[2, :, :],
              a   = bc.a[:, 1, :],
              b   = bc.b[:, 1, :],
              rhs =  rhs[:, 1, :], p=p, Δh = p.Δx)

    ## right
    @views apply_bc!(grid_ghost  = grid[nx,     :, :],
              grid_inside = grid[nx - 1, :, :],
              a   = bc.a[:, 2, :],
              b   = bc.b[:, 2, :],
              rhs =  rhs[:, 2, :], p=p, Δh = p.Δx)

    ## top
    @views apply_bc!(grid_ghost  = grid[:, ny,     :],
              grid_inside = grid[:, ny - 1, :],
              a   = bc.a[:, 3, :],
              b   = bc.b[:, 3, :],
              rhs =  rhs[:, 3, :], p=p, Δh = p.Δy)

    ## bottom
    @views apply_bc!(grid_ghost  = grid[:, 1, :],
              grid_inside = grid[:, 2, :],
              a   = bc.a[:, 4, :],
              b   = bc.b[:, 4, :],
              rhs =  rhs[:, 4, :], p=p, Δh = p.Δy)
end

function update_RHS_saturation!(; s_RHS::AbstractMatrix{T},
        P::Array{T, 1}, p::Parameters,
        t::Integer) where T<:AbstractFloat
    @views ρ̂₁_bound    = density.(ideal_gas,  P)
    @views ρ̂₂_bound    = density.(tait_C₅H₁₂, P)

    if t < p.relax_steps 
        ψ = (p.ψ - p.ψ₀) * t / relax_steps + p.ψ₀
    else
        ψ = p.ψ
    end

    @views s_RHS[:, 1] = saturation.(ρ̂₁_bound, ρ̂₂_bound, ψ,
                                         p.M[1], p.M[2])
    @views @. s_RHS[:, 2] = 1 - s_RHS[:, 1]
end

function update_RHS_saturation!(s_RHS::Array{T, 3}, P::Array{T, 2},
        p::Parameters, t::Integer) where T<:AbstractFloat
    ## bottom inlet
    @views update_RHS_saturation!(s_RHS = s_RHS[p.inlet, 4, :],
                                  P = fill(p.Pin, last(p.inlet) -
                                           first(p.inlet) + 1),
                           p = p, t = t)
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

function update_RHS_pressure!(P_RHS::Array{T, 3},
        p::Parameters, t::Integer) where T<:AbstractFloat
    if t < p.relax_steps
        P_inlet = (p.Pin - p.P₀) * t / p.relax_steps + p.P₀
    else
        P_inlet = p.Pin
    end

    ## bottom inlet
    @views @. P_RHS[p.inlet, 4, 1] = P_inlet

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

function apply_bc_pres_sat_den!(sys::System, p::Parameters,
        t::Integer)
    ## ∂P / ∂n = 0 at walls
    ## P = Pin     at inlet
    ## P = Pout    at outlet
    update_RHS_pressure!(sys.bc_RHS.P, p, t)
    P_tmp = zeros(p.nx, p.ny, 2)
    P_tmp[:, :, 1] = sys.P
    apply_bc!(P_tmp, p.bc.P, sys.bc_RHS.P, p)
    sys.P = P_tmp[:, :, 1]

    ## s_inlet = s_inlet(ρ̂₁(T, P), ρ̂₂(P, T), P, T)
    ## ∂s / ∂n = 0 everywhere else
    update_RHS_saturation!(sys.bc_RHS.s, sys.P, p, t)
    apply_bc!(sys.s, p.bc.s, sys.bc_RHS.s, p)

    ## ρ_k_bound = ρ̂_k_bound(T, P) * s_k_bound(T, P)
    update_RHS_density!(sys.bc_RHS.ρ, sys.s, sys.P, p)
    apply_bc!(sys.ρ, p.bc.ρ, sys.bc_RHS.ρ, p)
end

# function apply_bc_velocity!(sys::System, p::Parameters)
#     ## Darcy at inlet and outlet
#     ## v_n at walls
#     update_RHS_velocity!(sys.bc_RHS.v, sys.P, sys.s, p)
#     apply_bc!(sys.u, p.bc.u, sys.bc_RHS.u, p)
#     apply_bc!(sys.v, p.bc.v, sys.bc_RHS.v, p)
# end

function initial_conditions!(sys::System, p::Parameters)
    ## ∂P / ∂n = 0 at walls
    ## P = Pin     at inlet
    ## P = Pout    at outlet
    update_RHS_pressure!(sys.bc_RHS.P, p, 0)
    P_tmp = zeros(p.nx, p.ny, 2)
    P_tmp[:, :, 1] = sys.P
    apply_bc!(P_tmp, p.bc.P, sys.bc_RHS.P, p)
    sys.P = P_tmp[:, :, 1]

    initial_saturation!(sys, p)

    ## s_inlet = s_inlet(ρ̂₁(T, P), ρ̂₂(P, T), P, T)
    update_RHS_saturation!(sys.bc_RHS.s, sys.P, p, 0)
    apply_bc!(sys.s, p.bc.s, sys.bc_RHS.s, p)

    density!(sys)
    update_RHS_density!(sys.bc_RHS.ρ, sys.s, sys.P, p)
    apply_bc!(sys.ρ, p.bc.ρ, sys.bc_RHS.ρ, p)

    ## Darcy at inlet and outlet
    # update_RHS_velocity!(sys.bc_RHS.v, sys.P, sys.s, p)
    # apply_bc!(sys.u, p.bc.u, sys.bc_RHS.u, p)
    # apply_bc!(sys.v, p.bc.v, sys.bc_RHS.v, p)
end

function initial_RHS(p::Parameters)
    # u_RHS = zeros(max(p.nx, p.ny), 4, 2)
    # v_RHS = zeros(max(p.nx, p.ny), 4, 2)
    ρ_RHS = zeros(max(p.nx, p.ny), 4, 2)
    s_RHS = zeros(max(p.nx, p.ny), 4, 2)
    P_RHS = zeros(max(p.nx, p.ny), 4, 2)
    return BoundaryRHS(ρ_RHS, s_RHS, P_RHS)
end
# -------------------------------------------------------- #

# -------------------- Flux ------------------------------ #
function calculate_flux!(sys::System, p::Parameters)
    # ## Interpolation of densities
    # ## Φu = u_(i + 1/2, j) * (ρ_(i, j) + ρ_(i + 1, j)) / 2
    # for i = 1 : nx - 1
    #     @. @views sys.Φu[i, :, :] = sys.u[i, :, :] *
    #                                     (sys.ρ[i,     :, :] +
    #                                      sys.ρ[i + 1, :, :]) / 2
    # end

    # ## Φv = v_(i, j + 1/2) * (ρ_(i, j) + ρ_(i, j + 1)) / 2
    # for j = 1 : ny - 1
    #     @. @views sys.Φv[:, j, :] = sys.v[:, j, :] *
    #                                     (sys.ρ[:, j,     :] +
    #                                      sys.ρ[:, j + 1, :]) / 2
    # end

    ## Upwind scheme
    for i = 1 : nx
        for j = 1 : ny
            for k = 1 : 2
                if i < nx && sys.P[i, j] >= sys.P[i + 1, j]
                    ## ⟶ 
                    sys.Φu[i, j, k] = sys.u[i, j, k] *
                                      sys.ρ[i, j, k]
                elseif i < nx && sys.P[i, j] < sys.P[i + 1, j]
                    ## ⟵
                    sys.Φu[i, j, k] = sys.u[i, j, k] *
                                      sys.ρ[i + 1, j, k]
                end
                if j < ny && sys.P[i, j] >= sys.P[i, j + 1]
                    ## ↑
                    sys.Φv[i, j, k] = sys.v[i, j, k] *
                                      sys.ρ[i, j, k]
                elseif j < ny && sys.P[i, j] < sys.P[i, j + 1]
                    ## ↓
                    sys.Φv[i, j, k] = sys.v[i, j, k] *
                                      sys.ρ[i, j + 1, k]
                end
            end
        end
    end
end

## Gotta change all this!!
# function get_boundary(a::Array{T, 2},
#         j::Integer) where T<:AbstractFloat
#     @views a_bound = (a[:, j] + a[:, j + 1]) / 2
#     return a_bound
# end
# 
# function get_flux(ρ::Array{T, 2}, s::Array{T, 2}, v::Array{T, 2},
#         j::Integer) where T<:AbstractFloat
#     ρ_bound = get_boundary(ρ, j)
#     s_bound = get_boundary(s, j)
#     v_bound = get_boundary(v, j)
#     return ρ_bound ./ s_bound .* v_bound
# end
# 
# function get_total_flux(sys::System, p::Parameters)
#     flux_init  = zeros(2, p.nx)
#     flux_final = zeros(2, p.nx)
#     for k = 1 : 2
#         flux_init[k, :]  = get_flux(sys.ρ[:, :, k],
#                                     sys.s[:, :, k],
#                                     sys.v[:, :, k], 1)
#         flux_final[k, :] = get_flux(sys.ρ[:, :, k],
#                                     sys.s[:, :, k],
#                                     sys.v[:, :, k], p.ny - 1)
#     end
#     flux_init_total  = sum(flux_init,  dims = 2)
#     flux_final_total = sum(flux_final, dims = 2)
#     return flux_init_total, flux_final_total
# end

function get_content_inside(sys::System, p::Parameters)
    content = zeros(2)

    for k = 1 : 2
        content[k] = sum(sys.ρ[2 : nx - 1, 2 : ny - 1, k])
    end

    return content
end

function get_flux_in_out(sys::System, p::Parameters)
    flux_init  = zeros(2)
    flux_final = zeros(2)

    for k = 1 : 2
        flux_init[k]  = sum(sys.Φv[2 : nx - 1, 1, k])
        flux_final[k] = sum(sys.Φv[2 : nx - 1, ny, k])
    end

    return flux_init, flux_final
end
# -------------------------------------------------------- #
