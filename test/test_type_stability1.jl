## Testing for type stability
## This code is type stable
## Check: julia> @code_warntype test()

function f(g::Function, x::Int64)
     return g(x)
end

g1(x)::Function = x^2
g2(x)::Function = (1 - x)^2
const garr = (g1, g2)

function test()
    for k = 1 : 2
           f(garr[k], 20)
    end
end
