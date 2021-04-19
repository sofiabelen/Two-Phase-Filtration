include("structs.jl")

using PyPlot

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
    #scale = 100.0
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
