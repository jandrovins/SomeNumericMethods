#!/usr/bin/env julia
# Absolute error
function abs_error(x, y)
    return abs(x-y)
end

# Relative error, being x the newest value
function rel_error(x, y)
    return abs((x-y) / x)
end

function fixed_point(x₀::Float64, f::Function, g::Function, fₑ::Function, tolerance::Float64, iterations::Int64)
    xlist = Float64[x₀]
    ylist = Float64[f(x₀)]
    elist = Float64[Inf]

    # if x₀ is a root
    if f(x₀) == 0
        println("x₀ is a root!! x₀=$x₀")
        return x₀
    end
    
    x = x₀
    i = 0
    print("i  ,  Xᵢ  ,  f(Xᵢ)  ,  ")
    if fₑ == abs_error
        println(" Abs error")
    else
        println(" Rel error")
    end

    println("$i  ,  $(xlist[end])  ,  $(ylist[end])  ,  $(elist[end])")
    i+=1
    while elist[end]  > tolerance && i <= iterations
        x = g(x)
        append!(xlist, x)
        append!(ylist, f(x))
        append!(elist, fₑ(xlist[end], xlist[end-1]))
        println("$i  ,  $(xlist[end])  ,  $(ylist[end])  ,  $(elist[end])")
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
    return exp(x-10.5)+sin(x)+x^(3)-x-x
end

function gy(x)
    return cos(x)+exp(x-21/2)+3*x^2-2
end


x₀ = parse(Float64, ARGS[1])
tolerance = parse(Float64, ARGS[2])
iterations = parse(Int64, ARGS[3])
x = fixed_point(x₀, fy, gy, abs_error, tolerance, iterations)
