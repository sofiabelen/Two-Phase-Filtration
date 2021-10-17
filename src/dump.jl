include("structs.jl")

using PyPlot

function plot_density(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3},
        ρ::AbstractArray{T, 3},
        p::Parameters) where T<:AbstractFloat
    fig, axs = PyPlot.subplots(1, 2, figsize=(13, 5))
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
        axs[k].plot(rgx, v[:, y_index, k] .* ρ[:, y_index, k] .* s[:, y_index, k])
    end

    axs[1].set_title("Поток азота через границу y = " * string(y_value), y=1.08)
    axs[2].set_title("Поток пентана через границу y = " * string(y_value), y=1.08)
    
    savefig("../img/2phase-filtration-flux" * string(Integer(y_value)) * ".png", dpi=200)
    savefig("../img/2phase-filtration-flux" * string(Integer(y_value)) * ".svg")
end

function plot(sys::System, p::Parameters)
    plot_density(sys.u, sys.v, sys.ρ, p)
    plot_pressure(sys.u, sys.v, sys.P, p)
    plot_fluxes(sys.u, sys.v, sys.ρ, sys.s, p, 0.0, 1)
    plot_fluxes(sys.u, sys.v, sys.ρ, sys.s, p, 2.0, p.ny)
end

function dump(sys::System, p::Parameters)
    x = ideal_gas
    pressure1 = @. R * x.T * sys.ρ[:, :, 1] / p.M[1] / sys.s
    x = tait_C₅H₁₂
    ρ̂₂ = @. sys.ρ[:, :, 2] / (1 - sys.s)
    pressure2 = @. (x.B + x.P₀) * 10^((ρ̂₂ - x.ρ₀) / x.C / 
                                      ρ̂₂) - x.B

    writedlm("../dump/pressure.txt", sys.P[:, 1], ' ')
    writedlm("../dump/pressure1.txt", pressure1 .- sys.P, ' ')
    writedlm("../dump/pressure2.txt", pressure2 .- sys.P, ' ')
    writedlm("../dump/v1.txt", sys.v[2, :, 1], ' ')
    writedlm("../dump/v2.txt", sys.v[2, :, 2], ' ')
    writedlm("../dump/density1.txt", sys.ρ[2, :, 1], ' ')
    writedlm("../dump/density2.txt", sys.ρ[2, :, 2], ' ')
    writedlm("../dump/flow1.txt", sys.ρ[:, ny, 1] .*
             sys.s[:, ny, 1] .* sys.v[:, ny, 1], ' ')
    writedlm("../dump/flow2.txt", sys.ρ[:, ny, 2] .*
             sys.s[:, ny, 2] .* sys.v[:, ny, 2], ' ')
    writedlm("../dump/saturation.txt", sys.s[2, :, 1], ' ')
end
