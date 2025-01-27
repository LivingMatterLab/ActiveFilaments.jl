#######################################################################
### Mathematical structures relevant to the elephant trunk modeling ###
### expansion derived from the Active Filament theory.              ###
### See: https://doi.org/10.1101/2024.10.27.620533                  ###
#######################################################################
@enum Architecture longitudinal=1 helical=2 radial=3

struct PiecewiseFunction{T, P}
    objects::SVector{T, P}
    ranges::SVector{T, SVector{2, Float64}}
end

struct TrunkInterpolations{T}
    R1::SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}
    R2::SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}
    K::SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}
end

struct ClampingCondition
    r0::SVector{3, Float64}
    d10::SVector{3, Float64}
    d20::SVector{3, Float64}
    d30::SVector{3, Float64}
end

struct Sphere
    r::Float64
    c::SVector{3, Float64}
end

struct SphereJoint
    sphere::Sphere
    sD_hat::SVector{3, Float64}
    sdD_hat_1::SVector{3, Float64}
    sdD_hat_2::SVector{3, Float64}
    sdD_hat_3::SVector{3, Float64}
end

@with_kw struct Trunk{T, N}
    "Length of the trunk"
    L::Float64 = 1.0

    "Young's modulus"
    E::Float64
    "Poisson's ratio"
    ν::Float64

    "Lower material coordinate limits"
    Z1::SVector{T, Float64}
    "Upper material coordinate limits"
    Z2::SVector{T, Float64}

    "Material coordinate vectors"
    Z::SVector{T, LinRange{Float64, Int64}} = SVector{T, LinRange{Float64, Int64}}([LinRange(
                                                                                        Z1[i],
                                                                                        Z2[i],
                                                                                        N)
                                                                                    for i in 1:T])

    "Outer tapering angle (assuming linear tapering)"
    φ0::Float64

    "Outer radius of the filament at Z = 0"
    R00::Float64
    "Outer radius of the filament (at all Z in segment)"
    R0::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
        reduce(vcat, [transpose(R00 .- Z[i] * tan(φ0)) for i in 1:T]
    ),
    )

    "Inner radii array at Z = 0"
    R10::SMatrix{T, 5, Float64}
    "Inner radii array at Z = 0"
    R20::SMatrix{T, 5, Float64}
    "Inner radii array (at all Z in segment)"
    R1::SMatrix{T, 5, SVector{N, Float64}} = SMatrix{T, 5, SVector{N, Float64}}(
        [R10[i, j] .- R10[i, j] / R00 * Z[i] * tan(φ0) for i in 1:T, j in 1:5]
    )
    "Outer radii array (at all Z in segment)"
    R2::SMatrix{T, 5, SVector{N, Float64}} = SMatrix{T, 5, SVector{N, Float64}}(
        [R20[i, j] .- R20[i, j] / R00 * Z[i] * tan(φ0) for i in 1:T, j in 1:5]
    )

    "Outer tapering angle array (assuming linear tapering)"
    φ2::SMatrix{T, 5, Float64} = SMatrix{T, 5, Float64}(
        [atan(R20[i, j] / R00 * tan(φ0)) for i in 1:T, j in 1:5]
    )

    "Lower angular limits for architectures, right"
    Θ1R::SMatrix{T, 5, Float64}
    "Upper angular limits for architectures, right"
    Θ2R::SMatrix{T, 5, Float64}
    "Lower angular limits for architectures, left"
    Θ1L::SMatrix{T, 5, Float64} = 2 * pi .- Θ2R
    "Upper angular limits for architectures, left"
    Θ2L::SMatrix{T, 5, Float64} = 2 * pi .- Θ1R

    "Fiber architectures"
    architectures::SMatrix{T, 5, Architecture} = SMatrix{T, 5, Architecture}(
        reduce(vcat, [[longitudinal helical helical radial radial] for i in 1:T]
    ),
    )
    "Outer ventral helical angle (tilde)"
    α2_ovo::SVector{T, Float64}
    "Inner ventral helical angle (tilde)"
    α2_ivo::SVector{T, Float64}

    "Trunk stiffness K0"
    K0::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
        reduce(vcat, [transpose(E * pi * R0[i, :] .^ 2) for i in 1:T]),
    )
    "Trunk stiffness K1"
    K1::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
        reduce(vcat, [transpose(E * pi / 4.0 * R0[i, :] .^ 4) for i in 1:T]),
    )
    "Trunk stiffness K2"
    K2::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
        reduce(vcat, [transpose(E * pi / 4.0 * R0[i, :] .^ 4) for i in 1:T]),
    )
    "Trunk stiffness K3"
    K3::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
        reduce(vcat, [transpose(E * pi / (1.0 + ν) / 4.0 * R0[i, :] .^ 4) for i in 1:T]),
    )

    "Volumetric density"
    ρvol::Float64

    "Total mass (assumes linear tapering)"
    m::Float64 = ρvol * pi / 3.0 * L * (R00^2 + R00 * R0[end, end] + R0[end, end]^2)

    "Radial muscle contribution to bending"
    radial_muscle_bending::Bool = false

    "Default clamping location"
    r0::SVector{3, Float64} = @SVector [0.0, 0.0, 0.0]

    "Trunk 'joint' sphere"
    sphere_r::Float64
    xi_1::Float64 = 0.0
    xi_2::Float64 = pi / 2
    xi_d::Float64
    sphere_joint::SphereJoint = SphereJoint(
        Sphere(
            sphere_r,
            r0 -
            AngleAxis(xi_2, 0.0, 1.0, 0.0) * AngleAxis(xi_1, 1.0, 0.0, 0.0) *
            [0.0, 0.0, 1.0] * sphere_r
        ),
        AngleAxis(xi_2, 0.0, 1.0, 0.0) * AngleAxis(xi_1, 1.0, 0.0, 0.0) * [0.0, 0.0, 1.0],
        AngleAxis(xi_d, 0.0, 1.0, 0.0) * AngleAxis(xi_1, 1.0, 0.0, 0.0) * [1.0, 0.0, 0.0],
        AngleAxis(xi_d, 0.0, 1.0, 0.0) * AngleAxis(xi_1, 1.0, 0.0, 0.0) * [0.0, 1.0, 0.0],
        AngleAxis(xi_d, 0.0, 1.0, 0.0) * AngleAxis(xi_1, 1.0, 0.0, 0.0) * [0.0, 0.0, 1.0]
    )

    "Default clamping orientation"
    d10::SVector{3, Float64} = sphere_joint.sdD_hat_1
    d20::SVector{3, Float64} = sphere_joint.sdD_hat_2
    d30::SVector{3, Float64} = sphere_joint.sdD_hat_3
    clamping_condition::ClampingCondition = ClampingCondition(r0, d10, d20, d30)
    clamping_condition_unrolled::SVector{12, Float64} = [
        clamping_condition.r0...,
        clamping_condition.d10...,
        clamping_condition.d20...,
        clamping_condition.d30...
    ]
end

@with_kw struct TrunkFast{T, N}
    "Main trunk definition structure"
    trunk::Trunk{T, N}

    "Precomputed quantities"
    p::SVector{2, NTuple{5, SVector{4, SMatrix{T, N, Float64}}}} = compute_p(trunk)

    "Interpolating functions"
    interpolations::TrunkInterpolations = TrunkInterpolations(
        SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}(
            [PiecewiseFunction(
                 SVector{T, Interpolations.Extrapolation}(
                     [cubic_spline_interpolation(trunk.Z[i], trunk.R1[i, j]) for i in 1:T]
                 ),
                 SVector{T, SVector{2, Float64}}(
                     [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                 )
             )
             for j in 1:5]
        ),
        SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}(
            [PiecewiseFunction(
                 SVector{T, Interpolations.Extrapolation}(
                     [cubic_spline_interpolation(trunk.Z[i], trunk.R2[i, j]) for i in 1:T]
                 ),
                 SVector{T, SVector{2, Float64}}(
                     [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                 )
             )
             for j in 1:5]
        ),
        SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}(
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(
                         trunk.Z[i],
                         trunk.K1[i, :],
                         extrapolation_bc = Line()
                     ) for i in 1:T]
                ),
                SVector{T, SVector{2, Float64}}(
                    [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                )
            ),
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(
                         trunk.Z[i],
                         trunk.K2[i, :],
                         extrapolation_bc = Line()
                     ) for i in 1:T]
                ),
                SVector{T, SVector{2, Float64}}(
                    [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                )
            ),
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(
                         trunk.Z[i],
                         trunk.K3[i, :],
                         extrapolation_bc = Line()
                     ) for i in 1:T]
                ),
                SVector{T, SVector{2, Float64}}(
                    [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                )
            ),
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(
                         trunk.Z[i],
                         trunk.K0[i, :],
                         extrapolation_bc = Line()
                     ) for i in 1:T]
                ),
                SVector{T, SVector{2, Float64}}(
                    [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                )
            )
        )
    )

    "Density integral"
    L::Float64 = trunk.L
    R00::Float64 = trunk.R00
    R0_end::Float64 = trunk.R0[end, end]
    ρlin0Int::Function = Z -> begin
        pi / (3.0 * L^2) * trunk.ρvol *
        (
            (L - Z)^3 * R00^2 + (L - Z)^2 * (L + 2 * Z) * R00 * R0_end +
            (L^3 - Z^3) * R0_end^2
        )
    end
end

@with_kw struct ActivatedTrunkQuantities{T, N}
    trunkFast::TrunkFast{T, N}
    γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}

    u_hat_array = compute_uhat_array(trunkFast, γ)
    u_hat = compute_uhat_interpolations(trunkFast.trunk, u_hat_array)

    R_factor::Float64 = compute_R_factor(trunkFast.trunk, u_hat_array[4])
    new_K = SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}(
        PiecewiseFunction(
            SVector{T, Interpolations.Extrapolation}(
                [cubic_spline_interpolation(
                     trunkFast.trunk.Z[i],
                     R_factor^4 * trunkFast.trunk.K1[i, :],
                     extrapolation_bc = Line()
                 ) for i in 1:T]
            ),
            SVector{T, SVector{2, Float64}}(
                [SVector{2, Float64}([trunkFast.trunk.Z1[i], trunkFast.trunk.Z2[i]])
                 for
                 i in 1:T]
            )
        ),
        PiecewiseFunction(
            SVector{T, Interpolations.Extrapolation}(
                [cubic_spline_interpolation(
                     trunkFast.trunk.Z[i],
                     R_factor^4 * trunkFast.trunk.K2[i, :],
                     extrapolation_bc = Line()
                 ) for i in 1:T]
            ),
            SVector{T, SVector{2, Float64}}(
                [SVector{2, Float64}([trunkFast.trunk.Z1[i], trunkFast.trunk.Z2[i]])
                 for
                 i in 1:T]
            )
        ),
        PiecewiseFunction(
            SVector{T, Interpolations.Extrapolation}(
                [cubic_spline_interpolation(
                     trunkFast.trunk.Z[i],
                     R_factor^4 * trunkFast.trunk.K3[i, :],
                     extrapolation_bc = Line()
                 ) for i in 1:T]
            ),
            SVector{T, SVector{2, Float64}}(
                [SVector{2, Float64}([trunkFast.trunk.Z1[i], trunkFast.trunk.Z2[i]])
                 for
                 i in 1:T]
            )
        ),
        PiecewiseFunction(
            SVector{T, Interpolations.Extrapolation}(
                [cubic_spline_interpolation(
                     trunkFast.trunk.Z[i],
                     R_factor^2 * trunkFast.trunk.K0[i, :],
                     extrapolation_bc = Line()
                 ) for i in 1:T]
            ),
            SVector{T, SVector{2, Float64}}(
                [SVector{2, Float64}([trunkFast.trunk.Z1[i], trunkFast.trunk.Z2[i]])
                 for
                 i in 1:T]
            )
        )
    )
end