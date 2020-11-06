function divdif(x::Array{Float64,1},f::Array{Float64,1})
#
# Returns the vector of divided differences for the
# Newton form of the interpolating polynomial.
#
# ON ENTRY:
# x abscisses, given as a column vector;
# f ordinates, given as a column vector.
#
# ON RETURN:
# d divided differences
#
#
    n = length(x)
    d = deepcopy(f)
    for i=2:n
        for j=1:i-1
            d[i] = (d[j] - d[i])/(x[j] - x[i])
        end
    end
    return d
end



x = [1;1.2;1.4;1.6;1.8;2]
y = [0.6747;0.8491;1.1214;1.4921;1.9607;2.5258]
d = divdif(x, y)
println(d)
