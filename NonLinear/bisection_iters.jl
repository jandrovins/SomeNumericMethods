#!/usr/bin/env julia

bisection_iters(a, b, tol) = ceil( log2( (b-a)/tol )  )

a = parse(Float64, ARGS[1])
b = parse(Float64, ARGS[2])
tolerance = parse(Float64, ARGS[3])
iters = bisection_iters(a, b, tolerance)
println("Iterations: $iters")
