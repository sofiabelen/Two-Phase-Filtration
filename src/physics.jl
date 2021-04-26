include("structs.jl")

R = 8.314
MOLAR_MASS_IGAS = 0.028
MOLAR_MASS_PENTANE = 0.07215
tait_C₅H₁₂ = TaitEoS(35e6, 1e5, 0.2105, 616.18)
ideal_gas = IdealGasEoS(298.0, MOLAR_MASS_IGAS)
