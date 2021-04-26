include("structs.jl")

using PyPlot

function plot_density(ρ::AbstractMatrix{T}, p) where T<:AbstractFloat
    fig = PyPlot.figure(figsize=(9, 8))
    ax = PyPlot.axes()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    rgx = range(0, 2, length=p.nx)
    rgy = range(0, 2, length=p.ny)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))
    
    ## Velocity vector field
    #scale = 100.0
    quiver(x, y, u[:, :, 1]', v[:, :, 1]', color="r")
    quiver(x, y, u[:, :, 2]', v[:, :, 2]', color="b")
    
    pos = ax.contourf(X, Y, ρ',
        cmap=matplotlib.cm.viridis, alpha=0.5)
    fig.colorbar(pos, ax=ax)
    # cp = contour(X, Y, P', cmap=matplotlib.cm.viridis)

    PyPlot.title("2 Phase Filtration: Liquid and Ideal Gas")
    
    savefig("../img/2phase-filtration-density.png", dpi=200)
    # savefig("../img/2phase-filtration.svg")
end

function plot(u::AbstractArray{T, 3},
        v::AbstractArray{T, 3},
        P::AbstractMatrix{T}, p) where T<:AbstractFloat
    fig = PyPlot.figure(figsize=(9, 8))
    ax = PyPlot.axes()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
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
    
    ## Pressure contour map
    pos = ax.contourf(X, Y, P',
        cmap=matplotlib.cm.viridis, alpha=0.5)
    fig.colorbar(pos, ax=ax)
    # cp = contour(X, Y, P', cmap=matplotlib.cm.viridis)

    PyPlot.title("2 Phase Filtration: Liquid and Ideal Gas")
    
    savefig("../img/2phase-filtration.png", dpi=200)
    savefig("../img/2phase-filtration.svg")
end

function plot(sys::System, p::Parameters)
    plot(sys.u, sys.v, sys.P, p)
end

function dump(sys::System, p::Parameters)
    x = ideal_gas
    pressure1 = @. R * x.T * sys.ρ[:, :, 1] / p.M[1] / sys.s
    x = tait_C₅H₁₂
    ρ̂₂ = @. sys.ρ[:, :, 2] / (1 - sys.s)
    pressure2 = @. (x.B + x.P₀) * 10^((ρ̂₂ - x.ρ₀) / x.C / 
                                      ρ̂₂) - x.B

    writedlm("pressure.txt", sys.P[2, :], ' ')
    writedlm("pressure1.txt", pressure1 .- sys.P, ' ')
    writedlm("pressure2.txt", pressure2 .- sys.P, ' ')
    writedlm("v1.txt", sys.v[2, :, 1], ' ')
    writedlm("v2.txt", sys.v[2, :, 2], ' ')
    writedlm("density1.txt", sys.ρ[2, :, 1], ' ')
    writedlm("density2.txt", sys.ρ[2, :, 2], ' ')
    writedlm("flow1.txt", sys.ρ[2, :, 1] .*
             sys.s[2, :] .* sys.v[2, :, 1], ' ')
    writedlm("flow2.txt", sys.ρ[2, :, 2] .*
             (1 .- sys.s[2, :]) .* sys.v[2, :, 2], ' ')
    writedlm("saturation.txt", sys.s[2, :], ' ')
end
