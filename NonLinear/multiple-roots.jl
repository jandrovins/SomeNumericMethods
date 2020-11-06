#!/usr/bin/env julia
# Absolute error
function abs_error(x, y)
    return abs(x-y)
end

# Relative error, being x the newest value
function rel_error(x, y)
    return abs((x-y) / x)
end

function multiple_roots(x₀::Float64, f::Function, g::Function, h::Function,fₑ::Function, tolerance::Float64, iterations::Int64, α::Float64)
    xlist = Float64[x₀]
    ylist = Float64[f(x₀)]
    glist = Float64[g(x₀)]
    hlist = Float64[h(x₀)]
    elist = Float64[Inf]
    clist = Float64[Inf, Inf]

    # if x₀ is a root
    if f(x₀) == 0
        println("x₀ is a root!! x₀=$x₀")
        return x₀
    end
    
    x = x₀
    i = 0
    print("i , Xᵢ , f(Xᵢ) , g(Xᵢ) , h(Xᵢ)")
    if fₑ == abs_error
        print(", Abs error")
    else
        print(", Rel error")
    end
    println(", Convergence analysis")
    println("$i , $(xlist[end]) , $(ylist[end]) , $(glist[end]) , $(hlist[end]), $(elist[end]), $(clist[end])")
    i+=1
    while elist[end]  > tolerance && i <= iterations
        x = xlist[end] - ((ylist[end] * glist[end]) / (glist[end]^2 - ylist[end] * hlist[end]))
        append!(xlist, x)
        append!(ylist, f(x))
        append!(glist, g(x))
        append!(hlist, h(x))
        append!(elist, fₑ(xlist[end], xlist[end-1]))
        if i > 1
            append!(clist, elist[end] / elist[end-1]^α)
        end
        println("$i , $(xlist[end]) , $(ylist[end]) , $(glist[end]) , $(hlist[end]), $(elist[end]), $(clist[end])")
        i += 1
    end

    if elist[end] <= tolerance
        println("Converged!")
        return xlist[end]
    end
    println("Reached max iterations. Error: $(elist[end])")
    return xlist[end]
end

function fy(x)
    return (log(x) - 1)^3
end

function gy(x)
    return (3/x) * (log(x)-1)^2
end

function hy(x)
    return (-3 * log(x)^2 + 12 * log(x) - 9)/x^2
end

x₀ = parse(Float64, ARGS[1])
tolerance = parse(Float64, ARGS[2])
iterations = parse(Int64, ARGS[3])
α = parse(Float64, ARGS[4])
x = multiple_roots(x₀, fy, gy, hy, abs_error, tolerance, iterations, α)
