function f(x::Integer, y::Integer)
    return x, y
end

a = zeros(3, 3)
b = zeros(3, 3)

function g()
    for i = 1 : 3
        for j = 1 : 3
            a[i, j], b[i, j] = f(i, j)
        end
    end
end
println(a, b)
