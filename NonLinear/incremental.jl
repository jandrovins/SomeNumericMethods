#!/usr/bin/env julia
# Incremental method: returns x₁, x₂, success and xlist, where x₁ and x₂ contain a root.
# If the max number of iterations is surpassed, 'success' is false; else true. xlist is
# the list of all the values of x, including x₀. If a root is found, x₁ and x₂ will be the
# same value.
function incremental(x₀::Float64, Δx::Float64, iterations::Int64, f::Function)
    xlist = Float64[x₀, x₀]
    ylist = Float64[f(x₀), f(x₀)]
    if iterations <= 0 || typeof(iterations) != Int64
        println("Iterations must be and Integer64 greater than 0")
        return x₀, x₀, false, xlist, ylist
    end
    if f(x₀) == 0
        println("$x₀ is a root!")
        return x₀, x₀, true, xlist, ylist
    end

    x₁ = x₀
    i = 1
    println("Iteration , xᵢ , f(xᵢ)")
    println("0 , $(xlist[end]) , $(ylist[end])")
    while (ylist[end-1] * ylist[end] > 0) && (iterations != 0)
        # x₀ will be x₁, and x₁ will be (x₁ + Δx)
        x₀, x₁ = x₁, x₁ + Δx
        iterations -= 1
        append!(xlist, x₁)
        append!(ylist, f(x₁))
        println("$i , $(xlist[end]) , $(ylist[end])")
        i +=1
    end

    # if a root was found
    if f(x₁) == 0
        println("$x₁ is a root!")
        return x₁, x₁, true, xlist[2:end], ylist[2:end]
    end

    # if an interval for a root was found
    if ylist[end-1] * ylist[end] < 0
        println("A root was found between $x₀ and $x₁")
        return x₀, x₁, true, xlist[2:end], ylist[2:end]
    end

    # if nothing was successful
    return x₀, x₁, false, xlist[2:end], ylist[2:end]
end

function fy(x)
    return exp(-x)-sin(4*x)
end

start = parse(Float64, ARGS[1])
Δ = parse(Float64, ARGS[2])
iterations = parse(Int64, ARGS[3])
left, right, success, xlist, ylist =incremental(start, Δ, iterations, fy)
