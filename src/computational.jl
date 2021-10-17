# -------------------- Binary Search --------------------- #

function binary_search(; f::Function, a::T, b::T,
        niter=50, eps_x::T=1e-6) where T<:AbstractFloat
    fa = f(a)
    fb = f(b)
    mid = (a + b) / 2

    for i = 1 : niter
        mid = (a + b) / 2

        if abs(b - a) < eps_x
            return mid
        end

        fmid = f(mid)
        fmid == 0 && return mid

        if sign(fmid) == sign(fa)
            a = mid
            fa = fmid
        else
            b = mid
            fb = fmid
        end
    end
    return error("binary search fails to reach accuracy of ",
                 eps_x, " in ", niter, " iterations")
end
# -------------------------------------------------------- #

# -------------------- Secant Method --------------------- #

function secant_root_finder(; f::Function, a::T, b::T,
        niter=100, eps_x::T=1e-6) where T<:AbstractFloat

    if f(a) * f(b) >= 0
        return error("secant method fails:
                     condition f(a)f(b) < 0 not met.")
    end

    for i = 1 : niter
        x₀ = a - f(a) * (b - a) / (f(b) - f(a))
        f₀ = f(x₀)
        @debug x₀ f₀

        if abs(b - a) < eps_x
            return x₀
        end

        if f(a) * f₀ < 0
            b = x₀
        elseif f(b) * f₀ < 0
            a = x₀
        elseif f₀ == 0
            return x₀
        else
            return error("secant method fails.")
        end
    end

    return error("secant method fails to reach accuracy of ",
                 eps_x, " in ", niter, " iterations")
    # return a - f(a) * (b - a) / (f(b) - f(a))
end
# -------------------------------------------------------- #

# -------------------- Newton Raphson -------------------- #

cache_x₀ = 1e5
is_cache = 1

function newton_raphson(; f::Function, fder::Function,
    niter=10000, eps_x::T=1e-6)::T where T<:AbstractFloat

    x₀::AbstractFloat = 0
    x₁::AbstractFloat = 0
    ## TODO: what if not cached?
    if is_cache == 1
        x₀ = cache_x₀
    end

    xs = fill(0.0, niter)
    for i = 1 : niter
        @debug x₀ f(x₀)
        # println(x₀, " ", f(x₀))
        x₁ = x₀ - f(x₀) / fder(x₀)

        ## To avoid domain error: negative values in log
        if x₁ < eps_x
            x₁ = eps_x
        end

        if abs(x₁ - x₀) < eps_x || f(x₁) == 0
            is_cache = 1
            cache_x₀ = x₁
            return x₁
        end

        xs[i] = x₀
        x₀ = x₁
    end

    rgx = range(1e4, 1e7, length=10000)
    output = hcat(rgx, f.(rgx))
    writedlm("../dump/f.txt", output, ' ')
    writedlm("../dump/xs.txt", xs, ' ')

    return error("newton raphson method fails to reach accuracy of ",
                 eps_x, " in ", niter, " iterations Delta P = ",
                 abs(f(x₀)), "x₀ = ", x₀)
end
# -------------------------------------------------------- #
