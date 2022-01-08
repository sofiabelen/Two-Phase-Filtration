include("structs.jl")

using PyPlot

function extrapolate_velocities(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3}, p::Parameters) where T<:AbstractFloat
    u_extrap = zeros(nx - 2, ny - 2, 2)
    v_extrap = zeros(nx - 2, ny - 2, 2)

    for i = 1 : nx - 2
        @views @. u_extrap[i, :, :] = (u[i,     2 : ny - 1, :]  +
                                       u[i + 1, 2 : ny - 1, :]) / 2
    end
    
    for j = 1 : ny - 2
        @views @. v_extrap[:, j, :] = (v[2 : nx - 1, j, :] +
                                       v[2 : nx - 1, j + 1, :]) / 2
    end
    return u_extrap, v_extrap
end

function plot_density(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3},
        ρ::AbstractArray{T, 3},
        s::AbstractArray{T, 3},
        p::Parameters) where T<:AbstractFloat
    fig, axs = PyPlot.subplots(1, 2, figsize=(10, 4))
    axs[1].set_xlabel("x, m")
    axs[1].set_ylabel("y, m")

    rgx = range(0, 2, length = p.nx - 2)
    rgy = range(0, 2, length = p.ny - 2)
    x  = [rgx;]'
    y  = [rgy;]'
    X  = repeat([rgx;]', length(x))
    Y  = repeat([rgy;],1, length(y))
    
    axs[1].set_ylabel(L"$y, (m)$")

    @views ρ̂ = ρ ./s
    @views ρ̂_inside = ρ̂[2 : nx - 1, 2 : ny - 1, :]

    u_extrap, v_extrap = extrapolate_velocities(u, v, p)

    for k = 1 : 2
        axs[k].set_xlabel(L"$x, (m)$")
        axs[k].quiver(x, y, u_extrap[:, :, k]', v_extrap[:, :, k]')

        pos = axs[k].contourf(X, Y, ρ̂_inside[:, :, k]',
                              cmap=matplotlib.cm.viridis,
                              alpha=0.5)

        axs[k].contour(X, Y, ρ̂_inside[:, :, k]',
                       cmap=matplotlib.cm.viridis)

        fig.colorbar(pos, ax=axs[k],
                     label=L"$\rho, \left(\frac{kg}{m^3}\right)$")
    end

    # fig.suptitle("Профиль плотностей и поле скоростей")
    axs[1].set_title("Азот")
    axs[2].set_title("Пентан")
    
    savefig("../img/2phase-filtration-density.png", dpi=200)
    savefig("../img/2phase-filtration-density.svg")
end

function plot_pressure(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3},
        P::AbstractMatrix{T}, p) where T<:AbstractFloat
    fig = PyPlot.figure(figsize=(9, 8))
    ax = PyPlot.axes()
    
    rgx = range(0, 2, length=p.nx - 2)
    rgy = range(0, 2, length=p.ny - 2)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))

    u_extrap, v_extrap = extrapolate_velocities(u, v, p)
    P_inside = P[2 : nx - 1, 2 : ny - 1]
    
    ## Velocity vector field
    quiver(x, y, u_extrap[:, :, 1]', v_extrap[:, :, 1]', color="r")
    quiver(x, y, u_extrap[:, :, 2]', v_extrap[:, :, 2]', color="b")
    
    pos = ax.contourf(X, Y, P_inside',
        cmap=matplotlib.cm.viridis, alpha=0.5)
    fig.colorbar(pos, ax=ax, label="Pressure")
    
    ## Pressure contour map
    # cp = contour(X, Y, P', cmap=matplotlib.cm.viridis)

    PyPlot.title("2 Phase Filtration: Nitrogen and Pentane")
    
    savefig("../img/2phase-filtration.png", dpi=200)
    savefig("../img/2phase-filtration.svg")
end

# function plot_fluxes(u::AbstractArray{T, 3},
#         v::AbstractArray{T, 3},
#         ρ::AbstractArray{T, 3},
#         s::AbstractArray{T, 3},
#         p::Parameters,
#         y_value::T, y_index::Integer) where T<:AbstractFloat
#     fig, axs = PyPlot.subplots(1, 2, figsize=(13, 7))
#     rgx = range(0, 2, length=p.nx)
# 
#     for k = 1 : 2
#         axs[k].set_xlabel("x, (m)")
#         axs[k].set_ylabel(L"$\left( \rho \cdot v_y \right), \left( \frac{kg \cdot m}{s}\right)$")
#         flux = get_flux(ρ[:, :, k], s[:, :, k], v[:, :, k], y_index)
#         axs[k].plot(rgx, flux)
#     end
# 
#     axs[1].set_title("Поток азота через границу y = "   * string(y_value), y=1.08)
#     axs[2].set_title("Поток пентана через границу y = " * string(y_value), y=1.08)
#     
#     savefig("../img/2phase-filtration-flux" * string(Integer(y_value)) * ".png", dpi=200)
#     savefig("../img/2phase-filtration-flux" * string(Integer(y_value)) * ".svg")
# end

function plot(sys::System, p::Parameters)
    plot_density(sys.u, sys.v, sys.ρ, sys.s, p)
    plot_pressure(sys.u, sys.v, sys.P, p)
    # plot_fluxes(sys.u, sys.v, sys.ρ, sys.s, p, 0.0, 1)
    # plot_fluxes(sys.u, sys.v, sys.ρ, sys.s, p, 2.0, p.ny - 1)
end
