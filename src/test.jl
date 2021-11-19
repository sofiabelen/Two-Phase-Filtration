include("cfd.jl")

function check_pressure(sys::System)
    for index in CartesianIndices(sys.P)
        i, j = Tuple(index)
        if sys.P[i, j] < 0
            return error("Negative Pressure [i, j] = [", i, ", ",
                         j,"]")
        end
    end
end

function check_saturation(sys::System, p::Parameters)
    for k = 1 : 2
        for j = 2 : p.ny - 1
            for i = 2 : p.nx - 1
                if sys.s[i, j, k] < 0
                    return error("Negative saturation, component ",
                                 k, " [x, y] = [", i, ", ", j, "]")
                end
                if sys.s[i, j, k] > 1
                    return error("Saturation > 1, component ",
                                 k, " [x, y] = [", i, ", ", j, "]")
                end
            end
        end
    end
end

function check_density(sys::System, p::Parameters)
    for k = 1 : 2
        for j = 2 : p.ny - 1
            for i = 2 : p.nx - 1
                if sys.ρ[i, j, k] < 0
                    return error("Negative density, component ",
                                 k, " [x, y] = [", i, ", ", j, "]")
                end
            end
        end
    end
end

function check_bc(bc::BoundaryCondition, var::String)
    for k = 1 : 2
        for w = 1 : 4
            for x = 1 : min(nx, ny)
                if (bc.a[x, w, k]^2 + bc.b[x, w, k]^2) == 0
                    return error("Wrong ", var, " bc at k = ", k,
                                 " w = ", w, " x = ", x)
                end
            end
        end
    end
end

function check_bc(bc::BoundaryConditions{T}) where T<:AbstractFloat
    check_bc(bc.u, "u")
    check_bc(bc.v, "v")
    check_bc(bc.ρ, "ρ")
    check_bc(bc.s, "s")
    check_bc(bc.P, "P")
end

function check_continuity(syswork::System, sysnext::System,
        p::Parameters, t::Integer)
    flux_init_work, flux_final_work = get_total_flux(syswork, p)
    content_work = get_content_inside(syswork, p)
    flux_init_next, flux_final_next = get_total_flux(sysnext, p)
    content_next = get_content_inside(sysnext, p)

    influx_work  = flux_init_work  * p.Δx * p.Δt
    outflux_work = flux_final_work * p.Δx * p.Δt
    inside_work  = content_next    * p.Δx * p.Δy

    influx_next  = flux_init_next  * p.Δx * p.Δt
    outflux_next = flux_final_next * p.Δx * p.Δt
    inside_next  = content_next    * p.Δx * p.Δy

    mass_work = influx_work - outflux_work + inside_work
    mass_next = influx_next - outflux_next + inside_next
    # println(mass_work, mass_next)
    diff = abs.(mass_work - mass_next)

    if t == 1
        io = open("../dump/mass.txt", "w")
        write(io, "m1 m2\n")
    else
        io = open("../dump/mass.txt", "a+")
    end
    write(io, string(mass_work[1]), " ", string(mass_work[2]), "\n")

    # eps = 0.001
    # if any(x -> x > eps, abs.(mass_work - mass_next))
    #     return error("Continuity equation breaks ", diff)
    # end
end
