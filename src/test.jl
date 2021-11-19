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
    # check_bc(bc.u, "u")
    # check_bc(bc.v, "v")
    check_bc(bc.ρ, "ρ")
    check_bc(bc.s, "s")
    check_bc(bc.P, "P")
end

function check(syswork::System, sysnext::System,
        p::Parameters, t::Int)
        check_saturation(syswork, p)
        check_density(syswork, p)
        check_pressure(syswork)
end
