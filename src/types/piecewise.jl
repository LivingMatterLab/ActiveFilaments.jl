#########################################
### PiecewiseStructure implementation ###
#########################################
#region ===========================
# abstract type PiecewiseStructure end
"""
    $(TYPEDEF)

Piecewise structure used to define symbolic or numeric piecewise functions.

$(FIELDS)
"""
struct PiecewiseStructure{T} # Rename to PiecewiseFunction <: PiecewiseStructure after deserialization is no longer an issue
    "Expressions (symbolic or numeric) in each of the piecewise ranges"
    expressions::Vector{T}
    "Array of piecewise ranges (2-element vectors or tuples) for the argument"
    piecewiseRanges::Any
end

"""
    $(TYPEDSIGNATURES)

Addition operator overload for PiecewiseStructure + PiecewiseStructure.
"""
function Base.:+(x::PiecewiseStructure, y::PiecewiseStructure)
    expressions = []
    for (expression_x, expression_y) in zip(x.expressions, y.expressions)
        push!(expressions, expression_x + expression_y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Addition operator overload for Number + PiecewiseStructure.
"""
function Base.:+(x::Number, y::PiecewiseStructure)
    expressions = []
    for expression in y.expressions
        push!(expressions, x + expression)
    end
    PiecewiseStructure(expressions, y.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Addition operator overload for PiecewiseStructure + Number.
"""
function Base.:+(x::PiecewiseStructure, y::Number)
    expressions = []
    for expression in x.expressions
        push!(expressions, y + expression)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Negation operator overload for PiecewiseStructure.
"""
function Base.:-(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, -expression)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Subtraction operator overload for PiecewiseStructure - PiecewiseStructure.
"""
function Base.:-(x::PiecewiseStructure, y::PiecewiseStructure)
    expressions = []
    for (expression_x, expression_y) in zip(x.expressions, y.expressions)
        push!(expressions, expression_x - expression_y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Subtraction operator overload for Number - PiecewiseStructure.
The PiecewiseStructure is subtracted from number in each of
`y.piecewiseRanges` to produce another PiecewiseStructure.
"""
function Base.:-(x::Number, y::PiecewiseStructure)
    expressions = []
    for expression in y.expressions
        push!(expressions, x - expression)
    end
    PiecewiseStructure(expressions, y.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Subtraction operator overload for PiecewiseStructure - Number.
The Number is subtracted from the PiecewiseStructure in each of
`x.piecewiseRanges` to produce another PiecewiseStructure.
"""
function Base.:-(x::PiecewiseStructure, y::Number)
    expressions = []
    for expression in x.expressions
        push!(expressions, expression - y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Multiplication operator overload for PiecewiseStructure * PiecewiseStructure.
"""
function Base.:*(x::PiecewiseStructure, y::PiecewiseStructure)
    expressions = []
    for (expression_x, expression_y) in zip(x.expressions, y.expressions)
        push!(expressions, expression_x * expression_y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Multiplication operator overload for Number * PiecewiseStructure.
"""
function Base.:*(x::Number, y::PiecewiseStructure)
    expressions = []
    for expression in y.expressions
        push!(expressions, x * expression)
    end
    PiecewiseStructure(expressions, y.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Multiplication operator overload for PiecewiseStructure * Number.
"""
function Base.:*(x::PiecewiseStructure, y::Number)
    expressions = []
    for expression in x.expressions
        push!(expressions, expression * y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Division operator overload for PiecewiseStructure / PiecewiseStructure.
"""
function Base.:/(x::PiecewiseStructure, y::PiecewiseStructure)
    expressions = []
    for (expression_x, expression_y) in zip(x.expressions, y.expressions)
        push!(expressions, expression_x / expression_y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Division operator overload for Number / PiecewiseStructure.
"""
function Base.:/(x::Number, y::PiecewiseStructure)
    expressions = []
    for expression in y.expressions
        push!(expressions, x / expression)
    end
    PiecewiseStructure(expressions, y.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Division operator overload for PiecewiseStructure / Number.
"""
function Base.:/(x::PiecewiseStructure, y::Number)
    expressions = []
    for expression in x.expressions
        push!(expressions, expression / y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Power operator overload for PiecewiseStructure ^ PiecewiseStructure.
"""
function Base.:^(x::PiecewiseStructure, y::PiecewiseStructure)
    expressions = []
    for (expression_x, expression_y) in zip(x.expressions, y.expressions)
        push!(expressions, expression_x^expression_y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Power operator overload for Number ^ PiecewiseStructure.
"""
function Base.:^(x::Number, y::PiecewiseStructure)
    expressions = []
    for expression in y.expressions
        push!(expressions, x^expression)
    end
    PiecewiseStructure(expressions, y.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Power operator overload for PiecewiseStructure ^ Number.
"""
function Base.:^(x::PiecewiseStructure, y::Number)
    expressions = []
    for expression in x.expressions
        push!(expressions, expression^y)
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Log function overload for PiecewiseStructure.
"""
function Base.:log(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, log(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Sin function overload for PiecewiseStructure.
"""
function Base.:sin(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, sin(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Cos function overload for PiecewiseStructure.
"""
function Base.:cos(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, cos(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Tan function overload for PiecewiseStructure.
"""
function Base.:tan(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, tan(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Cot function overload for PiecewiseStructure.
"""
function Base.:cot(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, cot(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Sec function overload for PiecewiseStructure.
"""
function Base.:sec(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, sec(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Sqrt function overload for PiecewiseStructure.
"""
function Base.:sqrt(x::PiecewiseStructure)
    expressions = []
    for expression in x.expressions
        push!(expressions, sqrt(expression))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    $(TYPEDSIGNATURES)

Atan function overload for PiecewiseStructure.
"""
function Base.:atan(x::PiecewiseStructure, y::PiecewiseStructure)
    expressions = []
    for (expression_x, expression_y) in zip(x.expressions, y.expressions)
        push!(expressions, atan(expression_x, expression_y))
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end

"""
    (p::PiecewiseStructure)(x)

PiecewiseStructure functor. Evaluates the piecewise structure
at a given argument `x`.
"""
function (p::PiecewiseStructure)(x)
    for (expression, range) in zip(p.expressions, p.piecewiseRanges)
        if (x >= range[1] && x <= range[2])
            return isempty(methods(expression)) ? expression : expression(x)
        end
    end
end

"""
    $(FUNCTIONNAME)(x::PiecewiseStructure, arg; worldage = true)

Builds a PiecewiseStructure `x` with symbolic expressions as callable, 
fast functions with an argument `arg`.
"""
function eval_build_functions(x::PiecewiseStructure, arg; worldage = true)
    expressions = []
    for expression in x.expressions
        if worldage
            push!(
                expressions,
                eval(build_function(expression, arg, expression = Val{false}))
            )
        else
            push!(expressions, eval(build_function(expression, arg)))
        end
    end
    PiecewiseStructure(expressions, x.piecewiseRanges)
end
#endregion ===========================