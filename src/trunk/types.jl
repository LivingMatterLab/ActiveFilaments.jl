
@enum Architecture longitudinal = 1 helical = 2 radial = 3

struct PiecewiseFunction{T, P}
    objects::SVector{T, P}
    ranges::SVector{T, SVector{2, Float64}}
end

struct TrunkInterpolations{T}
    R1::SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}
    R2::SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}
    K::SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}
    # K0::PiecewiseFunction{T, Interpolations.Extrapolation}
    # K1::PiecewiseFunction{T, Interpolations.Extrapolation}
    # K2::PiecewiseFunction{T, Interpolations.Extrapolation}
    # K3::PiecewiseFunction{T, Interpolations.Extrapolation}
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
    # "Material coordinate vectors"
    # Z::SVector{T, SVector{N, Float64}} = 
    #     SVector{T, SVector{N, Float64}}([SVector{N, Float64}(LinRange(Z1[i], Z2[i], N)) for i in 1:T])
    
    "Material coordinate vectors"
    Z::SVector{T, LinRange{Float64, Int64}} = 
        SVector{T, LinRange{Float64, Int64}}([LinRange(Z1[i], Z2[i], N) for i in 1:T])
        
    "Outer tapering angle (assuming linear tapering)"
    φ0::Float64

    "Outer radius of the filament at Z = 0"
    R00::Float64
    "Outer radius of the filament (at all Z in segment)"
    R0::SMatrix{T, N, Float64} = 
        SMatrix{T, N, Float64}(
            reduce(vcat, [
                transpose(R00 .- Z[i] * tan(φ0)) for i in 1:T
            ]
            )
        )

    "Inner radii array at Z = 0"
    R10::SMatrix{T, 5, Float64}
    "Inner radii array at Z = 0"
    R20::SMatrix{T, 5, Float64}
    "Inner radii array (at all Z in segment)"
    R1::SMatrix{T, 5, SVector{N, Float64}} = 
        SMatrix{T, 5, SVector{N, Float64}}(
            [
                R10[i, j] .- R10[i, j] / R00 * Z[i] * tan(φ0) for i in 1:T, j in 1:5
            ]
        )
    "Outer radii array (at all Z in segment)"
    R2::SMatrix{T, 5, SVector{N, Float64}} = 
        SMatrix{T, 5, SVector{N, Float64}}(
            [
                R20[i, j] .- R20[i, j] / R00 * Z[i] * tan(φ0) for i in 1:T, j in 1:5
            ]
        )
    
    "Outer tapering angle array (assuming linear tapering)"
    φ2::SMatrix{T, 5, Float64} = 
        SMatrix{T, 5, Float64}(
            [
                atan(R20[i, j] / R00 * tan(φ0)) for i in 1:T, j in 1:5
            ]
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
    architectures::SMatrix{T, 5, Architecture} = 
        SMatrix{T, 5, Architecture}(
            reduce(vcat, [
                [longitudinal helical helical radial radial] for i in 1:T
            ]
            )
        )
    "Outer ventral helical angle (tilde)"
    α2_ovo::SVector{T, Float64}
    "Inner ventral helical angle (tilde)"
    α2_ivo::SVector{T, Float64}

    # "Helical angles (tilde)"
    # α2::SMatrix{T, 6, Union{Float64, Nothing}} = 
    #     SMatrix{T, 6, Union{Float64, Nothing}}(
    #         reduce(vcat, 
    #             [
    #                 [nothing, -α2_ovo, α2_ovo, α2_ivo, -α2_ivo, nothing] for i in 1:T
    #             ]
    #             )
    #     )
    
    "Trunk stiffness K0"
    K0::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
       reduce(vcat, [transpose(E * pi * R0[i, :] .^ 2) for i in 1:T])
    )
    "Trunk stiffness K1"
    K1::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
       reduce(vcat, [transpose(E * pi / 4.0 * R0[i, :] .^ 4) for i in 1:T])
    )
    "Trunk stiffness K2"
    K2::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
       reduce(vcat, [transpose(E * pi / 4.0 * R0[i, :] .^ 4) for i in 1:T])
    )
    "Trunk stiffness K3"
    K3::SMatrix{T, N, Float64} = SMatrix{T, N, Float64}(
       reduce(vcat, [transpose(E * pi / (1.0 + ν) / 4.0 * R0[i, :] .^ 4) for i in 1:T])
    )

    # "Precomputed quantities"
    # p::SVector{7, SVector{3, SMatrix{T, N, Float64}}} = compute_p(trunk)

    # "Interpolating functions"
    # interpolations::TrunkInterpolations = TrunkInterpolations(
    #     SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}(
    #         [PiecewiseFunction(
    #             SVector{T, Interpolations.Extrapolation}(
    #                 [cubic_spline_interpolation(Z[i], R1[i, j]) for i in 1:T]
    #             ),
    #             SVector{T, SVector{2, Float64}}(
    #                 [SVector{2, Float64}([Z1[i], Z2[i]]) for i in 1:T]
    #             )
    #         )
    #         for j in 1:5]
    #     ),
    #     SVector{5, PiecewiseFunction{T, Interpolations.Extrapolation}}(
    #         [PiecewiseFunction(
    #             SVector{T, Interpolations.Extrapolation}(
    #                 [cubic_spline_interpolation(Z[i], R2[i, j]) for i in 1:T]
    #             ),
    #             SVector{T, SVector{2, Float64}}(
    #                 [SVector{2, Float64}([Z1[i], Z2[i]]) for i in 1:T]
    #             )
    #         )
    #         for j in 1:5]
    #     ),
    #     PiecewiseFunction(
    #         SVector{T, Interpolations.Extrapolation}(
    #             [cubic_spline_interpolation(Z[i], K0[i, :]) for i in 1:T]
    #         ),
    #         SVector{T, SVector{2, Float64}}(
    #             [SVector{2, Float64}([Z1[i], Z2[i]]) for i in 1:T]
    #         )
    #     ),
    #     PiecewiseFunction(
    #         SVector{T, Interpolations.Extrapolation}(
    #             [cubic_spline_interpolation(Z[i], K1[i, :]) for i in 1:T]
    #         ),
    #         SVector{T, SVector{2, Float64}}(
    #             [SVector{2, Float64}([Z1[i], Z2[i]]) for i in 1:T]
    #         )
    #     ),
    #     PiecewiseFunction(
    #         SVector{T, Interpolations.Extrapolation}(
    #             [cubic_spline_interpolation(Z[i], K2[i, :]) for i in 1:T]
    #         ),
    #         SVector{T, SVector{2, Float64}}(
    #             [SVector{2, Float64}([Z1[i], Z2[i]]) for i in 1:T]
    #         )
    #     ),
    #     PiecewiseFunction(
    #         SVector{T, Interpolations.Extrapolation}(
    #             [cubic_spline_interpolation(Z[i], K3[i, :]) for i in 1:T]
    #         ),
    #         SVector{T, SVector{2, Float64}}(
    #             [SVector{2, Float64}([Z1[i], Z2[i]]) for i in 1:T]
    #         )
    #     )
    # )

    "Volumetric density"
    ρvol::Float64

    "Total mass (assumes linear tapering)"
    m::Float64 = ρvol * pi / 3.0 * L * (R00^2 + R00 * R0[end, end] + R0[end, end]^2)
end

@with_kw struct TrunkFast{T, N}
    "Main trunk definition structure"
    trunk::Trunk{T, N}

    "Precomputed quantities"
    p::SVector{2, NTuple{5, SVector{3, SMatrix{T, N, Float64}}}} = compute_p(trunk)

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
                [cubic_spline_interpolation(trunk.Z[i], trunk.K0[i, :], extrapolation_bc = Line()) for i in 1:T]
            ),
            SVector{T, SVector{2, Float64}}(
                [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
            )
            ),
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(trunk.Z[i], trunk.K1[i, :], extrapolation_bc = Line()) for i in 1:T]
                ),
                SVector{T, SVector{2, Float64}}(
                    [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                )
            ),
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(trunk.Z[i], trunk.K2[i, :], extrapolation_bc = Line()) for i in 1:T]
                ),
                SVector{T, SVector{2, Float64}}(
                    [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
                )
            ),
            PiecewiseFunction(
                SVector{T, Interpolations.Extrapolation}(
                    [cubic_spline_interpolation(trunk.Z[i], trunk.K3[i, :], extrapolation_bc = Line()) for i in 1:T]
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
    ρlin0Int::Function = Z -> 
        begin
            pi / (3.0 * L^2) * trunk.ρvol * ((L - Z)^3 * R00^2 + (L - Z)^2 * (L + 2 * Z) * R00 * R0_end + (L^3 - Z^3) * R0_end^2)
        end
end

@with_kw struct ActivatedTrunkQuantities{T, N}
    trunkFast::TrunkFast{T, N}
    γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}
    u_hat_array = compute_uhat_array(trunkFast, γ)
    u_hat = compute_uhat_interpolations(trunkFast.trunk, u_hat_array)
end