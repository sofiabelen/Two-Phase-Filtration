include("structs.jl")
include("cfd.jl")

function dump_total_liquid_flux(sys::System, p::Parameters,
        t::Integer)
    if t == 1
        io_init  = open("../dump/total_flux_init.txt",  "w")
        io_final = open("../dump/total_flux_final.txt", "w")
    else
        io_init  = open("../dump/total_flux_init.txt",  "a+")
        io_final = open("../dump/total_flux_final.txt", "a+")
    end

    writedlm(io_init,  sys.Φv[:, 1, 1] .\
                        sys.s[:, 1, 1], '\n')

    writedlm(io_final, sys.Φv[:, ny - 1, 2] .\
                        sys.s[:, ny - 1, 2], '\n')
end

function dump(sys::System, p::Parameters)
    # x = ideal_gas
    # pressure1 = @. R * x.T * sys.ρ[:, :, 1] / p.M[1] / sys.s
    # x = tait_C₅H₁₂
    # ρ̂₂ = @. sys.ρ[:, :, 2] / (1 - sys.s)
    # pressure2 = @. (x.B + x.P₀) * 10^((ρ̂₂ - x.ρ₀) / x.C / 
    #                                   ρ̂₂) - x.B

    @views ρ̂ = sys.ρ ./ sys.s

    writedlm("../dump/pressure.txt", sys.P[:, :], ' ')
    # writedlm("../dump/pressure1.txt", pressure1 .- sys.P, ' ')
    # writedlm("../dump/pressure2.txt", pressure2 .- sys.P, ' ')
    writedlm("../dump/u1.txt", sys.u[1 : nx - 1, :, 1], ' ')
    writedlm("../dump/u2.txt", sys.u[1 : nx - 1, :, 2], ' ')
    writedlm("../dump/v1.txt", sys.v[:, 1 : ny - 1, 1], ' ')
    writedlm("../dump/v2.txt", sys.v[:, 1 : ny - 1, 2], ' ')
    writedlm("../dump/density1.txt", ρ̂[:, :, 1], ' ')
    writedlm("../dump/density2.txt", ρ̂[:, :, 2], ' ')
    writedlm("../dump/flow1.txt", sys.Φv[:, :, 1] .\
                                  sys.s[:, :, 1], ' ')
    writedlm("../dump/flow2.txt", sys.Φv[:, :, 2] .\
                                  sys.s[:, :, 2], ' ')
    writedlm("../dump/saturation.txt", sys.s[:, :, 1], ' ')
end

function dump_mass(sys::System, p::Parameters, t::Integer)
    @views flux_init  = sum(sys.Φv[:, 1, :] ./
                            sys.s[:, 1, :], dims = 1)
    @views flux_final = sum(sys.Φv[:, ny - 1, :] ./
                            sys.s[:, ny - 1, :], dims = 1)

    content = get_content_inside(sys, p)
    flux_init  = reshape(flux_init,  size(content))
    flux_final = reshape(flux_final, size(content))

    influx  = flux_init  * p.Δx * p.Δt
    outflux = flux_final * p.Δx * p.Δt
    inside  = content    * p.Δx * p.Δy

    mass = influx - outflux + inside

    if t == 1
        io = open("../dump/mass.txt", "w")
        write(io, "m1 m2\n")
    else
        io = open("../dump/mass.txt", "a+")
    end

    write(io, string(mass[1]), " ", string(mass[2]), "\n")
end
