######### PiecewiseStructure integration methods
#region ===========================
"""
    $(TYPEDSIGNATURES)

Integrates the piecewise structure `x` with respect to `arg`.
"""
function integrate_piecewise(x::PiecewiseStructure, arg::Num)
    expressions = []
    for expression in x.expressions
        push!(expressions, integrate(expression, arg)[1])
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Integrates the piecewise structure with a symbolic antiderivative `ad` from `a` to `b`.
"""
function evaluate_integral_AD(ad::PiecewiseStructure, a::Number, b::Number)
    sign = 1.0
    if a > b
        temp = b
        b = a
        a = temp
        sign = -1.0
    end
    expressions_a = []
    ranges_a = []
    for (expression, range) in zip(ad.expressions, ad.piecewiseRanges)
        if (a >= range[1] && a <= range[2])
            if (b >= range[1] && b <= range[2])
                return sign * (expression(b) - expression(a))
            else
                push!(expressions_a, expression)
                push!(ranges_a, range)
            end
        elseif (b >= range[1] && b <= range[2])
            result = expression(b) - expression(range[1])
            for i in length(expressions_a):-1:2
                result += expressions_a[i](ranges_a[i][2]) -
                          expressions_a[i](ranges_a[i][1])
            end
            result += expressions_a[1](ranges_a[1][2]) - expressions_a[1](a)
            return sign * result
        end
    end
end
#endregion ===========================

######### Auxiliary integration methods
#region ===========================
"""
    $(TYPEDSIGNATURES)

Integrates the function with antiderivative `ad` from `a` to `b`.
"""
function evaluate_integral_AD(ad, a, b)
    sign = 1.0
    if a > b
        temp = b
        b = a
        a = temp
        sign = -1.0
    end

    return sign * (ad(b) - ad(a))
end
#endregion ===========================