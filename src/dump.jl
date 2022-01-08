include("structs.jl")
include("cfd.jl")

function dump_total_in_out_flux(sys::System, p::Parameters,
        t::Integer)
    if t == 1
        io_init  = open("../dump/total_flux_init.txt",  "w")
        io_final = open("../dump/total_flux_final.txt", "w")
        write(io_init,  "phi1 phi2\n")
        write(io_final, "phi1 phi2\n")
    else
        io_init  = open("../dump/total_flux_init.txt",  "a+")
        io_final = open("../dump/total_flux_final.txt", "a+")
    end

    influx, outflux = get_flux_in_out(sys, p) .* p.Δx .* p.Δt

    # fmt = Printf.Format("%.10g %.10g\n")
    # Printf.format(io_init,  fmt, influx[1],  influx[2])
    # Printf.format(io_final, fmt, outflux[2], outflux[2])
    write(io_init,  string(influx[1]), " ",  string(influx[2]), "\n")
    write(io_final, string(outflux[1])," ", string(outflux[2]), "\n")
    close(io_init)
    close(io_final)
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
    writedlm("../dump/flow1_in.txt",  sys.Φv[2 : nx - 1, 1,  1], ' ')
    writedlm("../dump/flow1_out.txt", sys.Φv[2 : nx - 1, ny - 1, 1], ' ')
    writedlm("../dump/flow2_in.txt",  sys.Φv[2 : nx - 1, 1,  2], ' ')
    writedlm("../dump/flow2_out.txt", sys.Φv[2 : nx - 1, ny - 1, 2], ' ')
    writedlm("../dump/flow1.txt", sys.Φv[:, :, 1], ' ')
    writedlm("../dump/flow2.txt", sys.Φv[:, :, 2], ' ')
    writedlm("../dump/saturation.txt", sys.s[:, :, 1], ' ')
end

function dump_mass(syswork::System, sysnext::System,
        p::Parameters, t::Integer)
    content_work = get_content_inside(syswork, p)
    content_next = get_content_inside(sysnext, p)

    flux_init, flux_final = get_flux_in_out(syswork, p)

    mass = zeros(2)
    influx = zeros(2)
    outflux = zeros(2)
    total_flux = zeros(2)
    side_flux1 = zeros(2)
    side_flux2 = zeros(2)

    for k = 1 : 2
        side_flux1[k] = sum(syswork.Φu[1,  2 : ny - 1, k]) * p.Δy * p.Δt
        side_flux2[k] = sum(syswork.Φu[nx, 2 : ny - 1, k]) * p.Δy * p.Δt
        total_flux[k] = sum(syswork.Φv[2 : nx - 1, :, k]) * p.Δx * p.Δt
        influx[k]  = flux_init[k]  * p.Δx * p.Δt
        outflux[k] = flux_final[k] * p.Δx * p.Δt
        inside  = (content_next[k] - content_work[k]) * p.φ * p.Δx * p.Δy
        # inside  = content_work[k] * p.Δx * p.Δy * p.φ
        mass[k] = inside - influx[k] + outflux[k]
    end

    if t == 1
        io = open("../dump/mass.txt", "w")
        write(io, "m1 m2 phi1 phi2 tot1 tot2 side_left1 side_right1 side_left2 side_right2\n")
    else
        io = open("../dump/mass.txt", "a+")
    end

    # write(io, string(mass[1]), " ", string(mass[2]), " ",
    #       string(influx[1]), " ", string(outflux[2])"\n")
    # Printf.Format("")
    fmt = Printf.Format("%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n")
    Printf.format(io, fmt, mass[1], mass[2], influx[1], influx[2],
                  total_flux[1], total_flux[2],
                  side_flux1[1], side_flux2[1],
                  side_flux1[2], side_flux2[2])
end
