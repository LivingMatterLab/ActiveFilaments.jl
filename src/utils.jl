mutable struct FiberID
    ringIndex::Number
    fiberIndex::Number
end

function angleBetweenVectors(u, v)
    return atan(norm(cross(u, v)), dot(u, v));
end

function projectOntoVector(u, v)
    return dot(u, v) / norm(v) ^ 2 * v;
end

function projectOntoPlane(u, n)
    return u - projectOntoVector(u, n);
end

function normalPlane3P(a, b, c)
    r1 = b - a;
    r2 = c - a;
    n = cross(r1, r2);
    return n / norm(n);
end

function rangeVectorAlwaysInclusive(start, stop, step)
    r = start:step:stop;
    if r[end] == stop
        return collect(r);   
    else
        return append!(collect(r), [stop]);
    end
end

function areaTriangle(x1, x2, x3)
    return 0.5 * norm(cross(x2 - x1, x3 - x1));
end

# https://en.wikipedia.org/wiki/Menger_curvature
function mengerCurvature(x1, x2, x3)
    A = areaTriangle(x1, x2, x3);
    denominator = norm(x1 - x2) * norm(x2 - x3) * norm(x3 - x1);
    return 4.0 * A / denominator;
end

function mengerRadius(x1, x2, x3)
    return 1.0 / mengerCurvature(x1, x2, x3);
end

function signedAngle(u, v, n)
    return atan(dot(cross(u, v), n), dot(u, v));
end