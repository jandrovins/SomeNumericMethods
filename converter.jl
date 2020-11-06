# Funciones!!
function times2(num)
    num_length = length(num)
    num = parse(BigInt, string(num))
    num_times2 = string(num * 2)
    num_times2_length = length(num_times2)
    leading_zeros = ""
    length_diff = num_length - num_times2_length
    if  length_diff > 0
        leading_zeros = "0"^length_diff
    end
    return leading_zeros * num_times2
end

function decimalstobinary(dec, places)
    dec = string(dec)
    result = ""
    for i in 1:places
        dec_length = length(dec)
        if parse(Int64, dec) == 0
            result *= "0"
            continue
        end

        dec_times2 = times2(dec)
        new_length = length(dec_times2)
        if new_length > dec_length
            result *= "1"
            dec = dec_times2[2:end]
        else
            result *= "0"
            dec = dec_times2
        end
    end
    return result
end

function ten_to(base, objective)
    intremainders = ""
    odigits = Dict(10 => "A", 11 => "B", 12 => "C", 13 => "D", 14=> "E", 15 => "F")
    sign = ""
    if base < 0
        sign = "-"
    end
    basedec, baseint = modf(abs(base))
    baseint = trunc(Int64, baseint)
    while baseint > 0
        remainder = baseint % objective
        if remainder >= 10 && haskey(odigits, remainder)
            remainder = odigits[remainder]
        else
            remainder = string(trunc(remainder))
        end
        intremainders = string(remainder) * intremainders
        baseint = trunc(Int64, baseint/objective)
    end
    return intremainders
end

function to_ten(digits, base)
    base10 = 0
    for i in string(digits)
        i = parse(Int64, i)
        base10 = base10 * base + i
    end
    return base10
end

function any_to_any(number, from, to)
    return ten_to(to_ten(number, from), to)
end

function float_dec_to_bin(num, size)
    num_split = split(string(num), ".")
    integer = ten_to(parse(Int64, num_split[1]), 2)
    integer_length = length(integer)
    size -= integer_length
    decimal = decimalstobinary(num_split[2], size)
    println("$num ⟹  $integer.$decimal")
    
    integer_length_bin = ten_to(integer_length, 2)
    println("$integer_length ⟹  $integer_length_bin")
end

float_dec_to_bin("0.972", 20)
#println(ten_to(488785, 2))
