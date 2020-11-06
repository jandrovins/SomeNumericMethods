#!/usr/bin/env julia
# Absolute error
function abs_error(x, y)
    return abs(x-y)
end

# Relative error, being x the newest value
function rel_error(x, y)
    return abs((x-y) / x)
end

# returns last x₀ and x₁ if a root was not found. If found, x₁ and x₀ will be equal
# 
function bisection(x₀::Float64, x₁::Float64, f::Function, fₑ::Function, tolerance::Float64, iterations::Int64)
# x₀ and x₁ are the first limits. f is the function to which we are finding  a root.
# fₑis the error funcion (abs or rel). iterations is the maximum number of iterations.
    # Assert that x₀ and x₁ have different signs
    if x₀ == x₁
        println("ERROR: el intervalo es el mismo número; $x₀ = $x₁")
    elseif f(x₀) * f(x₁) > 0 
        println("ERROR: los datos ingresados tienen el mismo signo")
        return -1
    elseif f(x₀) * f(x₁) ==0
        if f(x₀) == 0
            println("x₀ is the root!!")
        else
            println("x₁ is the root!!")
        end
    end

    leftlist = Float64[]
    midlist = Float64[]
    rightlist = Float64[]
    yleftlist = Float64[]
    ymidlist = Float64[]
    yrightlist = Float64[]
    elist = Float64[Inf]

    i = 0
    print("Iteration , x₀ , x₁ , xₘ , f(x₀) , f(x₁) , f(ₘ) ,")
    if fₑ == abs_error
        println(" Abs error")
    else
        println(" Rel error")
    end
    while f(x₁) != 0 && f(x₀) != 0 && i <= iterations && elist[end] > tolerance
        xₘ = (x₀ + x₁) / 2
        fx₀, fx₁,  fxₘ = f(x₀), f(x₁), f(xₘ)
        # Save x values
        append!(leftlist, x₀)
        append!(midlist, xₘ)
        append!(rightlist, x₁)
        # Save f(x) values
        append!(yleftlist, fx₀)
        append!(ymidlist, fxₘ)
        append!(yrightlist, fx₁)
        print("$(i+1) , $x₀ , $x₁ , $xₘ , $fx₀ , $fx₁ , $fxₘ")
        if fxₘ * fx₀ < 0
            x₁ = xₘ
        else
            x₀ = xₘ
        end

        if i != 0
            append!(elist, fₑ(midlist[end], midlist[end-1]))
        end

        i += 1
        println(" , $(elist[end])")
    end

    if f(x₀) * f(x₁) ==0
        if f(x₁) == 0
            println("The root is $x₁ !!")
            return x₁, x₁
        else
            println("The root is $x₀ !!")
            return x₀, x₀
        end
    end

    if elist[end] <= tolerance
        println("The tolerance $(tolerance) has been reached! x₀=$(midlist[end]), x₁=$(midlist[end]), xₘ=$(midlist[end]) and error=$(elist[end])")
        return midlist[end], midlist[end]
    end
    
    println("The root with tolerance $tolerance was not found. New limits: x₀=$x₀ , x₁=$x₁")
    return x₀, x₁
end

function fy(x)
    return log(x^2+x-0.6) + 0.5 -x
end

left = parse(Float64, ARGS[1])
right = parse(Float64, ARGS[2])
tolerance = parse(Float64, ARGS[3])
iterations = parse(Int64, ARGS[4])
left, right = bisection(left, right, fy, abs_error, tolerance, iterations)
