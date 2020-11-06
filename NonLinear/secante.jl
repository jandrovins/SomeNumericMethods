#!/usr/bin/env julia
# Absolute error
function abs_error(x, y)
    return abs(x-y)
end

# Relative error, being x the newest value
function rel_error(x, y)
    return abs((x-y) / x)
end

function fixed_point(x₀::Float64, x₁::Float64, f::Function, fₑ::Function, tolerance::Float64, iterations::Int64)
    xlist = Float64[x₀, x₁]
    ylist = Float64[f(x₀), f(x₁)]
    elist = Float64[Inf, Inf]

    # if x₀ is a root
    if f(x₀) == 0
        println("x₀ is a root!! x₀=$x₀")
        return x₀
    end
    # if x₁ is a root
    if f(x₁) == 0
        println("x₁ is a root!! x₁=$x₁")
        return x₁
    end
    
    print("i , Xᵢ , f(Xᵢ) , ")
    if fₑ == abs_error
        println(" Abs error")
    else
        println(" Rel error")
    end

    x = x₀
    i = 0
    println("$i , $(xlist[end-1]) , $(ylist[end-1]) , $(elist[end-1])")
    i+=1
    println("$i , $(xlist[end]) , $(ylist[end])  , $(elist[end])")
    i+=1
    while elist[end]  > tolerance && i <= iterations
        x = xlist[end] - ylist[end]*(xlist[end] - xlist[end-1])/(ylist[end]-ylist[end-1])
        append!(xlist, x)
        append!(ylist, f(x))
        append!(elist, fₑ(xlist[end], xlist[end-1]))
        println("$i , $(xlist[end]) , $(ylist[end])  , $(elist[end])")
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
    return exp(x) -2*x^8*exp(x/2)+x^16-0.1
end

x₀ = parse(Float64, ARGS[1])
x₁ = parse(Float64, ARGS[2])
tolerance = parse(Float64, ARGS[3])
iterations = parse(Int64, ARGS[4])
x = fixed_point(x₀, x₁, fy, abs_error, tolerance, iterations)
