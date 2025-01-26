# This is an estimate. This becomes exact with n -> Infinity
function estimateTotalBending(
    sol::ODESolution,
    filament,
    refVector;
    n = 20,
    normal = [1.0, 0.0, 0.0],
)
    Z = LinRange(0.0, filament.L, n)
    points = sol.(Z)

    totalBending = 0.0
    previous_point = points[1]
    for point in points
        totalBending += signedAngle(previous_point[10:12], point[10:12], normal)
        previous_point = point
    end

    return totalBending
end

function bendingAngle(
    sol::ODESolution,
    filament,
    refVector;
    largeAngle = false,
    fullRev = false,
    n_bending = 20,
    normal = [1.0, 0.0, 0.0],
)
    totalBending =
        estimateTotalBending(sol, filament, refVector, n = n_bending, normal = normal)
    tangent_end = sol[end][10:12]

    angle_end = signedAngle(refVector, tangent_end, normal)
    if totalBending >= 0
        if totalBending < pi
            return angle_end
        elseif totalBending < 2 * pi
            return pi + (pi + angle_end)
        else
            return 2 * pi + angle_end
        end
    else
        if totalBending > -pi
            return angle_end
        elseif totalBending > -2 * pi
            return -pi - (pi - angle_end)
        else
            return -2 * pi + angle_end
        end
    end
end

function bendingAngle(r::Vector, refVector)
    end_vector = [r[1][end] - r[1][end-1], r[2][end] - r[2][end-1], r[3][end] - r[3][end-1]]
    end_vector = end_vector / norm(end_vector)
    return end_vector[2] <= 0.0 ?
           atan(norm(cross(end_vector, refVector)), dot(end_vector, refVector)) :
           2 * pi - atan(norm(cross(end_vector, refVector)), dot(end_vector, refVector))
end

function loopRadius(sol::ODESolution, filament; f1 = 1.0 / 5.0, f2 = 2.0 / 5.0)
    point1 = sol(filament.L)[1:3]
    point2 = sol(filament.L - filament.L * f1)[1:3]
    point3 = sol(filament.L - filament.L * f2)[1:3]
    # Projection not needed for Menger curvature.

    radius = mengerRadius(point1, point2, point3)
    return radius
end

function loopRadius(r::Vector, L; f1 = 1.0 / 5.0, f2 = 2.0 / 5.0)
    n = size(r[1])[1]
    Z_range = LinRange(0.0, L, n)
    id1 = n
    id2 = findfirst(Z -> Z >= L * (1 - f1), Z_range)
    id3 = findfirst(Z -> Z >= L * (1 - f2), Z_range)

    rM = reduce(hcat, r)
    point1 = rM[id1, :]
    point2 = rM[id2, :]
    point3 = rM[id3, :]

    radius = mengerRadius(point1, point2, point3)
    return radius
end

function wrappingAngle(sol::ODESolution, refVector)
    end_vector = sol[end][10:12]
    end_vector_xy = [end_vector[1], end_vector[2], 0.0]
    if end_vector_xy == [0.0, 0.0, 0.0]
        return 0.0
    end
    end_vector_xy = end_vector_xy / norm(end_vector_xy)
    angle = atan(norm(cross(end_vector_xy, refVector)), dot(end_vector_xy, refVector))
    angle = (end_vector_xy[2] <= 0.0 ? angle : 2 * pi - angle)
    return angle
end

function tiltAngle(
    sol::ODESolution,
    filament,
    planeNormal;
    f1 = 1.0 / 5.0,
    f2 = 2.0 / 5.0,
    projectionPlaneNormal = [],
)
    point_end1 = sol(filament.L)[1:3]
    point_end2 = sol(filament.L - filament.L * f1)[1:3]
    point_end3 = sol(filament.L - filament.L * f2)[1:3]
    n = normalPlane3P(point_end2, point_end1, point_end3)

    if projectionPlaneNormal != []
        n = projectOntoPlane(n, projectionPlaneNormal)
        n = n / norm(n)
    end
    return angleBetweenVectors(n, planeNormal)
end

function tiltAngleCapped(
    sol::ODESolution,
    filament,
    planeNormal;
    f1 = 1.0 / 5.0,
    f2 = 2.0 / 5.0,
    kwargs...,
)
    angle = tiltAngle(sol, filament, planeNormal; f1 = f1, f2 = f2, kwargs...)
    if angle > pi / 2
        angle = tiltAngle(sol, filament, -planeNormal; f1 = f1, f2 = f2, kwargs...)
    end
    return angle
end

function tiltAngle(
    r::Vector,
    L,
    planeNormal;
    f1 = 1.0 / 5.0,
    f2 = 2.0 / 5.0,
    projectionPlaneNormal = [],
)
    n = size(r[1])[1]
    Z_range = LinRange(0.0, L, n)
    id1 = n
    id2 = findfirst(Z -> Z >= L * (1 - f1), Z_range)
    id3 = findfirst(Z -> Z >= L * (1 - f2), Z_range)

    rM = reduce(hcat, r)
    point_end1 = rM[id1, :]
    point_end2 = rM[id2, :]
    point_end3 = rM[id3, :]

    n = normalPlane3P(point_end2, point_end1, point_end3)
    println(n)
    if projectionPlaneNormal != []
        n = projectOntoPlane(n, projectionPlaneNormal)
        n = n / norm(n)
    end

    angle = angleBetweenVectors(n, planeNormal)
    if angle > pi / 2
        angle = tiltAngle(r, L, -planeNormal; f1 = f1, f2 = f2)
    end

    return angle
end

function fiberCurve(fiberID::FiberID, sol::ODESolution, Z, filament, activation_structure)
    θ0 = activation_structure[fiberID.ringIndex].θ0
    R2 = filament.rings[fiberID.ringIndex].geometry.R2
    α2 = filament.rings[fiberID.ringIndex].fiberArchitecture.α2
    arg = θ0 + Z * tan(α2) / R2
    return sol(Z)[1:3] + R2 * (cos(arg) * sol(Z)[4:6] + sin(arg) * sol(Z)[7:9])
end

### Limitations of fiberCurveD, fiberCurveRefD:
### - only for non-tapered filaments
### - currently works only if there is one fiber per ring. For this to work
###   with rings with multiple fibers, need to change θ0 to the angle for a
###   particular fiber in a ring given FiberID (not just θ0 for the whole ring)
### These limitations propagate to: fiberLength, fiberLengthRef, fiberStrain
function fiberCurveD(fiberID::FiberID, sol::ODESolution, Z, filament, activation_structure)
    θ0 = activation_structure[fiberID.ringIndex].θ0
    R2 = filament.rings[fiberID.ringIndex].geometry.R2
    α2 = filament.rings[fiberID.ringIndex].fiberArchitecture.α2
    if (θ0 isa PiecewiseStructure)
        θ0 = θ0(Z)
    end
    if (α2 isa PiecewiseStructure)
        α2 = α2(Z)
    end
    tana = tan(α2)
    arg = θ0 + Z * tana / R2
    c = cos(arg)
    s = sin(arg)
    solVal = sol(Z)
    d1 = solVal[4:6]
    d2 = solVal[7:9]

    solD = sol(Z, Val{1})
    rD = solD[1:3]
    d1D = solD[4:6]
    d2D = solD[7:9]

    return c * d2 * tana + R2 * c * d1D + s * (-d1 * tana + R2 * d2D) + rD
end

function fiberCurveRefD(fiberID::FiberID, Z, filament, activation_structure)
    θ0 = activation_structure[fiberID.ringIndex].θ0
    R2 = filament.rings[fiberID.ringIndex].geometry.R2
    α2 = filament.rings[fiberID.ringIndex].fiberArchitecture.α2
    if (θ0 isa PiecewiseStructure)
        θ0 = θ0(Z)
    end
    if (α2 isa PiecewiseStructure)
        α2 = α2(Z)
    end
    tana = tan(α2)
    arg = θ0 + Z * tana / R2
    c = cos(arg)
    s = sin(arg)
    d1 = [1.0, 0.0, 0.0]
    d2 = [0.0, 1.0, 0.0]
    rD = [0.0, 0.0, 1.0]

    return c * d2 * tana + s * (-d1 * tana) + rD
end

function fiberLength(fiberID::FiberID, sol::ODESolution, filament, activation_structure)
    integral, error = quadgk(
        Z -> norm(fiberCurveD(fiberID, sol, Z, filament, activation_structure)),
        0.0, filament.L)
    return integral
end

function fiberLengthRef(fiberID::FiberID, filament, activation_structure)
    integral, error = quadgk(
        Z -> norm(fiberCurveRefD(fiberID, Z, filament, activation_structure)),
        0.0, filament.L)
    return integral
end

function fiberStrain(fiberID::FiberID, sol::ODESolution, filament, activation_structure)
    return fiberLength(fiberID, sol, filament, activation_structure) /
           fiberLengthRef(fiberID, filament, activation_structure) - 1.0
end

