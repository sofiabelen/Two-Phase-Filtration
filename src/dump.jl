include("structs.jl")

using PyPlot

function get_boundary(a::Array{T, 2},
        i::Integer) where T<:AbstractFloat
    a_bound::Array{AbstractFloat} = zeros(size(a))
    @. @views a_bound = (a[:, i] + a[:, i + 1]) / 2
    return a_bound
end

function get_flux(ρ::Array{T, 2}, s::Array{T, 2}, v::Array{T, 2},
        i::Integer) where T<:AbstractFloat
    ρ_bound = get_boundary(ρ[:, :], i)
    s_bound = get_boundary(s[:, :], i)
    v_bound = get_boundary(v[:, :], i)
    return ρ_bound .* s_bound .* v_bound
end

function plot_density(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3},
        ρ::AbstractArray{T, 3},
        p::Parameters) where T<:AbstractFloat
    fig, axs = PyPlot.subplots(1, 2, figsize=(10, 4))
    axs[1].set_xlabel("x, m")
    axs[1].set_ylabel("y, m")

    rgx = range(0, 2, length=p.nx)
    rgy = range(0, 2, length=p.ny)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))
    
    axs[1].set_ylabel(L"$y, (m)$")
    for k = 1 : 2
        axs[k].set_xlabel(L"$x, (m)$")
        axs[k].quiver(x, y, u[:, :, k]', v[:, :, k]')
        pos = axs[k].contourf(X, Y, ρ[:, :, k]',
                              cmap=matplotlib.cm.viridis,
                              alpha=0.5)
        axs[k].contour(X, Y, ρ[:, :, k]',
                              cmap=matplotlib.cm.viridis)
        fig.colorbar(pos, ax=axs[k], label=L"$\rho, \left(\frac{kg}{m^3}\right)$")
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
    
    rgx = range(0, 2, length=p.nx)
    rgy = range(0, 2, length=p.ny)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))
    
    ## Velocity vector field
    # scale = 100.0
    quiver(x, y, u[:, :, 1]', v[:, :, 1]', color="r")
    quiver(x, y, u[:, :, 2]', v[:, :, 2]', color="b")
    
    pos = ax.contourf(X, Y, P',
        cmap=matplotlib.cm.viridis, alpha=0.5)
    fig.colorbar(pos, ax=ax, label="Pressure")
    
    ## Pressure contour map
    # cp = contour(X, Y, P', cmap=matplotlib.cm.viridis)

    PyPlot.title("2 Phase Filtration: Nitrogen and Pentane")
    
    savefig("../img/2phase-filtration.png", dpi=200)
    savefig("../img/2phase-filtration.svg")
end

function plot_fluxes(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3},
        ρ::AbstractArray{T, 3},
        s::AbstractArray{T, 3},
        p::Parameters,
        y_value::T, y_index::Integer) where T<:AbstractFloat
    fig, axs = PyPlot.subplots(1, 2, figsize=(13, 7))
    rgx = range(0, 2, length=p.nx)

    for k = 1 : 2
        axs[k].set_xlabel("x, (m)")
        axs[k].set_ylabel(L"$\left( \rho \cdot v_y \right), \left( \frac{kg \cdot m}{s}\right)$")
        flux = get_flux(ρ[:, :, k], s[:, :, k], v[:, :, k], y_index)
        # not needed anymore:
        # axs[k].plot(rgx, v[:, y_index, k] .* ρ[:, y_index, k] .* s[:, y_index, k])
        axs[k].plot(rgx, flux)
    end

    axs[1].set_title("Поток азота через границу y = "   * string(y_value), y=1.08)
    axs[2].set_title("Поток пентана через границу y = " * string(y_value), y=1.08)
    
    savefig("../img/2phase-filtration-flux" * string(Integer(y_value)) * ".png", dpi=200)
    savefig("../img/2phase-filtration-flux" * string(Integer(y_value)) * ".svg")
end

function dump_total_liquid_flux(sys::System, p::Parameters,
        t::Integer)
    # total_flux_final = 0
    # total_flux_init  = 0

    # ρ₂_init  = get_boundary(sys.ρ[:, :, 2], 1)
    # ρ₂_final = get_boundary(sys.ρ[:, :, 2], p.nx - 1)
    # s₂_init  = get_boundary(sys.s[:, :, 2], 1)
    # s₂_final = get_boundary(sys.s[:, :, 2], p.nx - 1)
    # v₂_init  = get_boundary(sys.v[:, :, 2], 1)
    # v₂_final = get_boundary(sys.v[:, :, 2], p.nx - 1)
    flux_init  = get_flux(sys.ρ[:, :, 2], sys.s[:, :, 2],
                         sys.v[:, :, 2], 1)
    flux_final = get_flux(sys.ρ[:, :, 2], sys.s[:, :, 2],
                         sys.v[:, :, 2], p.ny - 1)

    total_flux_init  = sum(flux_init)
    total_flux_final = sum(flux_final)

    if t == 1
        io_init  = open("../dump/total_flux_init.txt",  "w")
        io_final = open("../dump/total_flux_final.txt", "w")
    else
        io_init  = open("../dump/total_flux_init.txt",  "a+")
        io_final = open("../dump/total_flux_final.txt", "a+")
    end
    writedlm(io_init,  total_flux_init,  '\n')
    writedlm(io_final, total_flux_final, '\n')
end

function plot(sys::System, p::Parameters)
    plot_density(sys.u, sys.v, sys.ρ, p)
    plot_pressure(sys.u, sys.v, sys.P, p)
    # plot_fluxes(sys.u, sys.v, sys.ρ, sys.s, p, 0.0, 1)
    # plot_fluxes(sys.u, sys.v, sys.ρ, sys.s, p, 2.0, p.ny - 1)
end

function dump(sys::System, p::Parameters)
    # x = ideal_gas
    # pressure1 = @. R * x.T * sys.ρ[:, :, 1] / p.M[1] / sys.s
    # x = tait_C₅H₁₂
    # ρ̂₂ = @. sys.ρ[:, :, 2] / (1 - sys.s)
    # pressure2 = @. (x.B + x.P₀) * 10^((ρ̂₂ - x.ρ₀) / x.C / 
    #                                   ρ̂₂) - x.B

    writedlm("../dump/pressure.txt", sys.P[:, :], ' ')
    # writedlm("../dump/pressure1.txt", pressure1 .- sys.P, ' ')
    # writedlm("../dump/pressure2.txt", pressure2 .- sys.P, ' ')
    writedlm("../dump/v1.txt", sys.v[:, :, 1], ' ')
    writedlm("../dump/v2.txt", sys.v[:, :, 2], ' ')
    writedlm("../dump/density1.txt", sys.ρ[:, :, 1], ' ')
    writedlm("../dump/density2.txt", sys.ρ[:, :, 2], ' ')
    writedlm("../dump/flow1.txt", sys.ρ[:, :, 1] .*
             sys.s[:, :, 1] .* sys.v[:, :, 1], ' ')
    writedlm("../dump/flow2.txt", sys.ρ[:, :, 2] .*
             sys.s[:, :, 2] .* sys.v[:, :, 2], ' ')
    writedlm("../dump/saturation.txt", sys.s[:, :, 1], ' ')
end
