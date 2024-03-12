function (p::PiecewiseFunction{T, Float64} where T)(x::Float64)
    for (object, range) in zip(p.objects, p.ranges)
        if (x >= range[1] && x <= range[2])
            return object;
        end
    end
end

function (p::PiecewiseFunction{T, Interpolations.Extrapolation} where T)(x::Float64)
    for (object, range) in zip(p.objects, p.ranges)
        if (x >= range[1] && x <= range[2])
            return object(x);
        end
    end
end

