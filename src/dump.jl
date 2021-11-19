include("structs.jl")
include("cfd.jl")

function dump_total_liquid_flux(sys::System, p::Parameters,
        t::Integer)
    flux_init  = get_flux(sys.ρ[:, :, 2],
                          sys.s[:, :, 2],
                          sys.v[:, :, 2], 1)
    flux_final = get_flux(sys.ρ[:, :, 2],
                          sys.s[:, :, 2],
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
    writedlm("../dump/v1.txt", sys.v[:, :, 1], ' ')
    writedlm("../dump/v2.txt", sys.v[:, :, 2], ' ')
    writedlm("../dump/density1.txt", ρ̂[:, :, 1], ' ')
    writedlm("../dump/density2.txt", ρ̂[:, :, 2], ' ')
    writedlm("../dump/flow1.txt", sys.ρ[:, :, 1] .\
             sys.s[:, :, 1] .* sys.v[:, :, 1], ' ')
    writedlm("../dump/flow2.txt", sys.ρ[:, :, 2] .\
             sys.s[:, :, 2] .* sys.v[:, :, 2], ' ')
    writedlm("../dump/saturation.txt", sys.s[:, :, 1], ' ')
end
