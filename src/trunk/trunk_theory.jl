@enum Side right = 1 left = 2

function (p::PiecewiseFunction{T, Float64} where T)(x::Float64)
    for (object, range) in zip(p.objects, p.ranges)
        if (x >= range[1] && x <= range[2])
            return object;
        end
    end
end

function (p::PiecewiseFunction{T, Interpolations.Extrapolation} where T)(x)
    for (object, range) in zip(p.objects, p.ranges)
        if (x >= range[1] && x <= range[2]) # Effectively: [a, b] for first segment, and then (a, b] for every subsequent
            return object(x);
        end
    end

    if (x > p.ranges[end][2])
        return p.objects[end](x)
    else
        return p.objects[1](x)
    end
end

# For plotting
function (p::PiecewiseFunction{T, Interpolations.Extrapolation} where T)(i::Integer)
    return p.objects[i]
end

function constant_activation(trunk::Trunk{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}) where {T, N}
    γR = γ[1]
    γL = γ[2]

    AR_1 = MMatrix{T, 5, Float64}(undef)
    AL_1 = MMatrix{T, 5, Float64}(undef)
    AR_2 = MMatrix{T, 5, Float64}(undef)
    AL_2 = MMatrix{T, 5, Float64}(undef)
    AR_3 = MMatrix{T, 5, Float64}(undef)
    AL_3 = MMatrix{T, 5, Float64}(undef)
    AR_0 = MMatrix{T, 5, Float64}(undef)
    AL_0 = MMatrix{T, 5, Float64}(undef)
    for i in 1:T
        for j in 1:5
            Θ1R = trunk.Θ1R[i, j]
            Θ2R = trunk.Θ2R[i, j]
            Θ1L = trunk.Θ1L[i, j]
            Θ2L = trunk.Θ2L[i, j]

            AR_1[i, j] = γR[i, j] * (cos(Θ1R) - cos(Θ2R))
            AL_1[i, j] = γL[i, j] * (cos(Θ1L) - cos(Θ2L))

            AR_2[i, j] = γR[i, j] * (-sin(Θ1R) + sin(Θ2R))
            AL_2[i, j] = γL[i, j] * (-sin(Θ1L) + sin(Θ2L))

            AR_3[i, j] = γR[i, j] * (-Θ1R + Θ2R)
            AL_3[i, j] = γL[i, j] * (-Θ1L + Θ2L)

            AR_0[i, j] = γR[i, j] * (-Θ1R + Θ2R)
            AL_0[i, j] = γL[i, j] * (-Θ1L + Θ2L)
        end
    end
    AR = (AR_1, AR_2, AR_3, AR_0)
    AL = (AL_1, AL_2, AL_3, AL_0)

    AR, AL
end

# These deltas assume constant linear tapering
function deltas_longitudinal(trunk::Trunk{T, N}, Z_index::Integer, muscle_index::Integer) where {T, N}
    ν = trunk.ν

    R1 = trunk.R1[Z_index, muscle_index]
    R2 = trunk.R2[Z_index, muscle_index]
    c_phi = tan(trunk.φ2[Z_index, muscle_index]) ./ R2

    δ1 = ν / 3.0 * (R1.^3 - R2.^3) + (1 + ν) ./ (c_phi.^2) .* (R2 - R1) + (1 + ν) ./ (c_phi.^3) .* (atan.(c_phi .* R1) - atan.(c_phi .* R2))
    δ2 = δ1
    δ3 = SVector{N, Float64}(zeros(N))

    δ0 = 0.5 * (
        (R1.^2 - R2.^2) * ν -
        (1 + ν) * log.(
            (1 .+ c_phi.^2 .* R1.^2) ./ (1 .+ c_phi.^2 .* R2.^2)
            ) ./ c_phi.^2
        )

    δ1, δ2, δ3, δ0
end

function deltas_helical(trunk::Trunk{T, N}, Z_index::Integer, muscle_index::Integer, side::Side) where {T, N}
    ν = trunk.ν

    R1 = trunk.R1[Z_index, muscle_index]
    R2 = trunk.R2[Z_index, muscle_index]

    c_phi = tan(trunk.φ2[Z_index, muscle_index]) ./ R2

    if muscle_index == 2
        α2 = (side == right ? -trunk.α2_ovo[Z_index] : trunk.α2_ovo[Z_index])
    else muscle_index 
        α2 = (side == right ? trunk.α2_ivo[Z_index] : -trunk.α2_ivo[Z_index])
    end
    c_a = tan(α2) ./ R2

    B22 = sqrt.(1 .+ c_phi.^2 .* R2.^2)
    B21 = sqrt.(1 .+ c_phi.^2 .* R1.^2)
    A = sqrt.(Complex.(c_phi.^2 - c_a.^2))

    δ1 = (R1.^3 - R2.^3) * ν / 3.0 + (1 + ν) ./ (c_a .* c_phi .* (c_a.^2 - c_phi.^2)) .*
            (c_phi .* atan.(c_a .* R1) - c_a .* atan.(c_phi .* R1) 
                - c_phi .* atan.(c_a .* R2) + c_a .* atan.(c_phi .* R2))
    δ2 = δ1
    δ3 = (B22 - B21) ./ (c_a .* c_phi.^2) + (atan.(c_a .* B21 ./ A) - atan.(c_a .* B22 ./ A)) ./ (c_a.^2 .* A)

    δ0 = 0.5 * (
        (R1.^2 - R2.^2) * ν + 
        (1 + ν) * log.(
            ((1 .+ c_phi.^2 .* R1.^2) .* (1 .+ c_a.^2 .* R2.^2)) ./
            ((1 .+ c_a.^2 .* R1.^2) .* (1 .+ c_phi.^2 .* R2.^2))
            ) ./
            (c_a.^2 - c_phi.^2)
        )

    δ1, δ2, real(δ3), δ0
end

function deltas_radial(trunk::Trunk{T, N}, Z_index::Integer, muscle_index::Integer) where {T, N}
    ν = trunk.ν

    R1 = trunk.R1[Z_index, muscle_index]
    R2 = trunk.R2[Z_index, muscle_index]

    δ1 = trunk.radial_muscle_bending ? ν / 3.0 * (R1.^3 - R2.^3) : SVector{N, Float64}(zeros(N))
    δ2 = δ1
    δ3 = SVector{N, Float64}(zeros(N))

    δ0 = 0.5 * (R1.^2 - R2.^2) * ν

    δ1, δ2, δ3, δ0
end

function compute_p(trunk::Trunk{T, N}) where {T, N}
    E = trunk.E
    K1 = trunk.K1
    K2 = trunk.K2
    K3 = trunk.K3
    K = (trunk.K1, trunk.K2, trunk.K3, trunk.K0)
    
    # p_dl_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # p_ovo_R_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # p_ovo_L_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # p_ivo_R_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # p_ivo_L_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # p_dr_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # p_vr_V = MVector{T, MVector{3, SVector{N, Float64}}}
    # for Z_index in 1:T
    #     δ_dl = deltas_longitudinal(trunk, Z_index, 1)
    #     δ_ovo_R = deltas_helical(trunk, Z_index, 2, right)
    #     δ_ovo_L = deltas_helical(trunk, Z_index, 2, left)
    #     δ_ivo_R = deltas_helical(trunk, Z_index, 3, right)
    #     δ_ivo_L = deltas_helical(trunk, Z_index, 3, left)
    #     δ_dr = deltas_radial(trunk, Z_index, 4)
    #     δ_vr = deltas_radial(trunk, Z_index, 5)

    #     p_dl_V[Z_index] = E * @MVector [δ_dl[i] ./ K[i][Z_index, :] for i in 1:3]
    #     p_ovo_R_V[Z_index] = E * @MVector [δ_ovo_R[i] ./ K[i][Z_index, :] for i in 1:3]
    #     p_ovo_L_V[Z_index] = E * @MVector [δ_ovo_L[i] ./ K[i][Z_index, :] for i in 1:3]
    #     p_ivo_R_V[Z_index] = E * @MVector [δ_ivo_R[i] ./ K[i][Z_index, :] for i in 1:3]
    #     p_ivo_L_V[Z_index] = E * @MVector [δ_ivo_L[i] ./ K[i][Z_index, :] for i in 1:3]
    #     p_dr_V[Z_index] = E * @MVector [δ_dr[i] ./ K[i][Z_index, :] for i in 1:3]
    #     p_vr_V[Z_index] = E * @MVector [δ_vr[i] ./ K[i][Z_index, :] for i in 1:3]

    #     # for i in 1:3
    #     #         E * 
    #     #         δ_lo[i] +
    #     #         δ_ovo_R[i] * AR[i][Z_index, 2] +
    #     #         δ_ovo_L[i] * AL[i][Z_index, 2] +  
    #     #         δ_ivo_R[i] * AR[i][Z_index, 3] +
    #     #         δ_ivo_L[i] * AL[i][Z_index, 3] +
    #     #         δ_dr[i] * (AR[i][Z_index, 4] + AL[i][Z_index, 4]) +
    #     #         δ_vr[i] * (AR[i][Z_index, 5] + AL[i][Z_index, 5])
    #     # end
    # end

    # # Reshape
    # p_dl = MVector{3, MMatrix{T, N, Float64}}
    # p_ovo_R = MVector{3, MMatrix{T, N, Float64}}
    # p_ovo_L = MVector{3, MMatrix{T, N, Float64}}
    # p_ivo_R = MVector{3, MMatrix{T, N, Float64}}
    # p_ivo_L = MVector{3, MMatrix{T, N, Float64}}
    # p_dr = MVector{3, MMatrix{T, N, Float64}}
    # p_vr = MVector{3, MMatrix{T, N, Float64}}
    # for i in 1:3
    #     for Z_index in 1:T
    #         p_dl[i][Z_index, :] = p_dl_V[Z_index][i]
    #         p_ovo_R[i][Z_index, :] = p_ovo_R_V[Z_index][i]
    #         p_ovo_L[i][Z_index, :] = p_ovo_L_V[Z_index][i]
    #         p_ivo_R[i][Z_index, :] = p_ivo_R_V[Z_index][i]
    #         p_ivo_L[i][Z_index, :] = p_ivo_L_V[Z_index][i]
    #         p_dr[i][Z_index, :] = p_dr_V[Z_index][i]
    #         p_vr[i][Z_index, :] = p_vr_V[Z_index][i]
    #     end
    # end
    
    δ_dl = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    δ_ovo_R = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    δ_ovo_L = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    δ_ivo_R = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    δ_ivo_L = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    δ_dr = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    δ_vr = Vector{NTuple{4, SVector{N, Float64}}}(undef, T)
    for Z_index in 1:T
        δ_dl[Z_index] = deltas_longitudinal(trunk, Z_index, 1)
        δ_ovo_R[Z_index] = deltas_helical(trunk, Z_index, 2, right)
        δ_ovo_L[Z_index] = deltas_helical(trunk, Z_index, 2, left)
        δ_ivo_R[Z_index] = deltas_helical(trunk, Z_index, 3, right)
        δ_ivo_L[Z_index] = deltas_helical(trunk, Z_index, 3, left)
        δ_dr[Z_index] = deltas_radial(trunk, Z_index, 4)
        δ_vr[Z_index] = deltas_radial(trunk, Z_index, 5)
    end

    # p_dl = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ovo_R = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ovo_L = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ivo_R = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ivo_L = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_dr = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_vr = MVector{3, MMatrix{T, N, Float64}}(undef)

    
    # for i in 1:3
    #     for Z_index in 1:T
    #         println(typeof(E * δ_dl[Z_index][i] ./ K[i][Z_index, :]))
    #         println(typeof(p_dl[i][Z_index, :]))
    #         p_dl[i][Z_index, :] = E * δ_dl[Z_index][i] ./ K[i][Z_index, :]
    #         p_ovo_R[i][Z_index, :] = E * δ_ovo_R[Z_index][i] ./ K[i][Z_index, :]
    #         p_ovo_L[i][Z_index, :] = E * δ_ovo_L[Z_index][i] ./ K[i][Z_index, :]
    #         p_ivo_R[i][Z_index, :] = E * δ_ivo_R[Z_index][i] ./ K[i][Z_index, :]
    #         p_ivo_L[i][Z_index, :] = E * δ_ivo_L[Z_index][i] ./ K[i][Z_index, :]
    #         p_dr[i][Z_index, :] = E * δ_dr[Z_index][i] ./ K[i][Z_index, :]
    #         p_vr[i][Z_index, :] = E * δ_vr[Z_index][i] ./ K[i][Z_index, :]
    #     end
    # end

    # println([reduce(hcat, E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T for i in 1:3])
    # println(typeof([reduce(hcat, E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T for i in 1:3]))
    # println(size([[reduce(vcat, E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T] for i in 1:3][1]))
    # println([[reduce(vcat, E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T] for i in 1:3][1])
    # println(typeof([[reduce(hcat, E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T] for i in 1:3][1]))

    p_dl = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    p_ovo_R = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_ovo_R[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    p_ovo_L = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_ovo_L[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    p_ivo_R = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_ivo_R[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    p_ivo_L = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_ivo_L[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    p_dr = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_dr[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    p_vr = @SVector [SMatrix{T, N, Float64}(reduce(vcat, [transpose(E * δ_vr[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T])) for i in 1:4]
    
    # p_dl = MVector{3, MMatrix{T, N, Float64}}([reduce(hcat, E * δ_dl[Z_index][i] ./ K[i][Z_index, :]) for Z_index in 1:T] )
    # p_ovo_R = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ovo_L = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ivo_R = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_ivo_L = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_dr = MVector{3, MMatrix{T, N, Float64}}(undef)
    # p_vr = MVector{3, MMatrix{T, N, Float64}}(undef)

    @SVector [(p_dl, p_ovo_R, p_ivo_R, p_dr, p_vr), (p_dl, p_ovo_L, p_ivo_L, p_dr, p_vr)]
    # @SVector [(p_dl, p_dl), (p_ovo_R, p_ovo_L), (p_ivo_R, p_ivo_L), (p_dr, p_dr), (p_vr, p_vr)]
    # @SVector [p_dl, p_ovo_R, p_ovo_L, p_ivo_R, p_ivo_L, p_dr, p_vr]
end

function compute_uhat_array(trunkFast::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}})  where {T, N}
    A = constant_activation(trunkFast.trunk, γ)
    p = trunkFast.p;
    u_hat = @SVector [
                        SMatrix{T, N, Float64}(
                        reduce(vcat, 
                            [
                            transpose(
                            p[1][1][u][i, :] * A[1][u][i, 1] + p[2][1][u][i, :] * A[2][u][i, 1] +
                            p[1][2][u][i, :] * A[1][u][i, 2] + p[2][2][u][i, :] * A[2][u][i, 2] +
                            p[1][3][u][i, :] * A[1][u][i, 3] + p[2][3][u][i, :] * A[2][u][i, 3] +
                            p[1][4][u][i, :] * A[1][u][i, 4] + p[2][4][u][i, :] * A[2][u][i, 4] +
                            p[1][5][u][i, :] * A[1][u][i, 5] + p[2][5][u][i, :] * A[2][u][i, 5]
                            )
                            for i in 1:T
                            ]
                            )
                        ) * (u == 2 ? -1.0 : 1.0) .+ (u == 4 ? 1.0 : 0.0)
                    for u in 1:4
                    ]

    u_hat
end

function compute_uhat_interpolations(trunk::Trunk{T, N}, u_hat_array::SArray) where {T, N}
    u_hat = 
        SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}(
            [
            PiecewiseFunction(
            SVector{T, Interpolations.Extrapolation}(
                [cubic_spline_interpolation(trunk.Z[i], u_hat_array[u][i, :], extrapolation_bc = Line()) for i in 1:T]
            ),
            SVector{T, SVector{2, Float64}}(
                [SVector{2, Float64}([trunk.Z1[i], trunk.Z2[i]]) for i in 1:T]
            )
            )
            for u in 1:4]
        )
    
    u_hat
end

# Assumes incompressibility and constant ζ_hat in each segment
function compute_R_factor(trunk::Trunk{T, N}, ζ_hat_array::SMatrix{T, N, Float64}) where {T, N}
    R_factor = sqrt(trunk.L / sum([ζ_hat_array[i, 1] * (trunk.Z2[i] - trunk.Z1[i]) for i in 1:T]))
    return R_factor
end

# R_factor in the current configuration (with external loading)
function compute_R_factor_current(trunk::Trunk{T, N}, sol; n = 200) where {T, N}
    Z_for_len = LinRange(0.0, trunk.L, n)
    len = 0.0
    r_Z = [sol(Z_for_len[i])[1:3] for i in 1:n]
    for i in 1:(n - 1)
        len = len + euclidean(r_Z[i], r_Z[i + 1])
    end
    R_factor = sqrt(trunk.L / len)
    println(R_factor)
    return R_factor
end

function intrinsic_trunk_de_SA(u, p, Z)
    u1_hat = p[1](Z)
    u2_hat = p[2](Z)
    u3_hat = p[3](Z)
    ζ_hat = p[4](Z)
    # ζ_hat = 1.0
    
    du1 = ζ_hat * u[10];
    du2 = ζ_hat * u[11];
    du3 = ζ_hat * u[12];
    du4 = ζ_hat * (u3_hat * u[7] - u2_hat * u[10]);
    du5 = ζ_hat * (u3_hat * u[8] - u2_hat * u[11]);
    du6 = ζ_hat * (u3_hat * u[9] - u2_hat * u[12]);
    du7 = ζ_hat * (u1_hat * u[10] - u3_hat * u[4]);
    du8 = ζ_hat * (u1_hat * u[11] - u3_hat * u[5]);
    du9 = ζ_hat * (u1_hat * u[12] - u3_hat * u[6]);
    du10 = ζ_hat * (u2_hat * u[4] - u1_hat * u[7]);
    du11 = ζ_hat * (u2_hat * u[5] - u1_hat * u[8]);
    du12 = ζ_hat * (u2_hat * u[6] - u1_hat * u[9]);
    SVector{12}(du1, du2, du3, du4, du5, du6, du7, du8, du9, du10, du11, du12)
end

function solveIntrinsic(trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}},
                        u0 = SVector{12, Float64}([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]), 
                        Zspan = (0.0, trunk.trunk.L); 
                        kwargs...) where {T, N}
    activation = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
    p = activation.u_hat
    println(typeof(p))
    
    prob = ODEProblem(intrinsic_trunk_de_SA, u0, Zspan, p);
    
    sol = solve(prob, AutoVern7(Rodas4()), dt = trunk.trunk.L / 100.0, abstol = 1e-12, reltol = 1e-12);

    sol, activation
end


function self_weight_trunk_de!(du, u, 
        p::Tuple{Float64, Function, SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}, SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}, SVector{3, Float64}, SVector{12, Float64}}, Z) where T
    u1_hat = p[4][1](Z)
    u2_hat = p[4][2](Z)
    u3_hat = p[4][3](Z)
    ζ_hat = p[4][4](Z)
    # ζ_hat = 1.0
    
    ρlinInt = p[2](Z);
    
    F1 = dot(p[5], u[4:6])
    F2 = dot(p[5], u[7:9])
    F3 = dot(p[5], u[10:12])
    n1 = -p[1] * ρlinInt * u[6] + F1;
    n2 = -p[1] * ρlinInt * u[9] + F2;
    n3 = -p[1] * ρlinInt * u[12] + F3;

    ζF = n3 / p[3][4](Z) + 1.0;
    ζ = ζ_hat * ζF;
    u1 = u1_hat + u[13] / p[3][1](Z);
    u2 = u2_hat + u[14] / p[3][2](Z);
    u3 = u3_hat + u[15] / p[3][3](Z);

    du[1] = ζ * u[10];
    du[2] = ζ * u[11];
    du[3] = ζ * u[12];
    du[4] = ζ_hat * (u3 * u[7] - u2 * u[10]);
    du[5] = ζ_hat * (u3 * u[8] - u2 * u[11]);
    du[6] = ζ_hat * (u3 * u[9] - u2 * u[12]);
    du[7] = ζ_hat * (u1 * u[10] - u3 * u[4]);
    du[8] = ζ_hat * (u1 * u[11] - u3 * u[5]);
    du[9] = ζ_hat * (u1 * u[12] - u3 * u[6]);
    du[10] = ζ_hat * (u2 * u[4] - u1 * u[7]);
    du[11] = ζ_hat * (u2 * u[5] - u1 * u[8]);
    du[12] = ζ_hat * (u2 * u[6] - u1 * u[9]);
    du[13] = ζ_hat * (u3 * u[14] - u2 * u[15]) + ζ * n2;
    du[14] = ζ_hat * (u1 * u[15] - u3 * u[13]) - ζ * n1;
    du[15] = ζ_hat * (u2 * u[13] - u1 * u[14]);
end

function self_weight_wrap_nbodyc_trunk_de!(du, u, 
        p::Tuple{
            Float64, 
            Function, 
            SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}, 
            SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}, 
            SVector{3, Float64}, 
            SVector{12, Float64}
            }, Z) where T
    u1_hat = p[4][1](Z)
    u2_hat = p[4][2](Z)
    u3_hat = p[4][3](Z)
    ζ_hat = p[4][4](Z)
    
    w = p[5][1]
    L = p[5][2]
    Z_0 = p[5][3]
    w_n = w * (L - max(Z, Z_0))

    ρlinInt = p[2](Z);
    n_factor = -p[1] * ρlinInt + w_n

    n1 = n_factor * u[6];
    n2 = n_factor * u[9];
    n3 = n_factor * u[12];

    ζF = n3 / p[3][4](Z) + 1.0;
    ζ = ζ_hat * ζF;
    u1 = u1_hat + u[13] / p[3][1](Z);
    u2 = u2_hat + u[14] / p[3][2](Z);
    u3 = u3_hat + u[15] / p[3][3](Z);

    du[1] = ζ * u[10];
    du[2] = ζ * u[11];
    du[3] = ζ * u[12];
    du[4] = ζ_hat * (u3 * u[7] - u2 * u[10]);
    du[5] = ζ_hat * (u3 * u[8] - u2 * u[11]);
    du[6] = ζ_hat * (u3 * u[9] - u2 * u[12]);
    du[7] = ζ_hat * (u1 * u[10] - u3 * u[4]);
    du[8] = ζ_hat * (u1 * u[11] - u3 * u[5]);
    du[9] = ζ_hat * (u1 * u[12] - u3 * u[6]);
    du[10] = ζ_hat * (u2 * u[4] - u1 * u[7]);
    du[11] = ζ_hat * (u2 * u[5] - u1 * u[8]);
    du[12] = ζ_hat * (u2 * u[6] - u1 * u[9]);
    du[13] = ζ_hat * (u3 * u[14] - u2 * u[15]) + ζ * n2;
    du[14] = ζ_hat * (u1 * u[15] - u3 * u[13]) - ζ * n1;
    du[15] = ζ_hat * (u2 * u[13] - u1 * u[14]);
end

function self_weight_wrap_trunk_de!(du, u, 
        p::Tuple{
            Float64, 
            Function, 
            SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}, 
            SVector{4, PiecewiseFunction{T, Interpolations.Extrapolation}}, 
            Tuple{Float64, Float64, Float64, SVector{3, Float64}, Function}, 
            SVector{12, Float64}
            }, Z) where T
    u1_hat = p[4][1](Z)
    u2_hat = p[4][2](Z)
    u3_hat = p[4][3](Z)
    ζ_hat = p[4][4](Z)
    
    w = p[5][1]
    L = p[5][2]
    Z_0 = p[5][3]
    r_C = p[5][4]
    R0_act = p[5][5](Z)

    w_n = w * (L - max(Z, Z_0))

    moment_r = r_C - u[1:3]
    moment_r = moment_r / norm(moment_r) * R0_act
    l1 = w * (u[4] * moment_r[2] - u[5] * moment_r[1])
    l2 = w * (u[7] * moment_r[2] - u[8] * moment_r[1])
    l3 = w * (u[10] * moment_r[2] - u[11] * moment_r[1])

    # l1 = 0
    # l2 = 0
    # l3 = 0

    ρlinInt = p[2](Z);
    n_factor = -p[1] * ρlinInt + w_n

    n1 = n_factor * u[6];
    n2 = n_factor * u[9];
    n3 = n_factor * u[12];

    ζF = n3 / p[3][4](Z) + 1.0;
    ζ = ζ_hat * ζF;
    u1 = u1_hat + u[13] / p[3][1](Z);
    u2 = u2_hat + u[14] / p[3][2](Z);
    u3 = u3_hat + u[15] / p[3][3](Z);

    du[1] = ζ * u[10];
    du[2] = ζ * u[11];
    du[3] = ζ * u[12];
    du[4] = ζ_hat * (u3 * u[7] - u2 * u[10]);
    du[5] = ζ_hat * (u3 * u[8] - u2 * u[11]);
    du[6] = ζ_hat * (u3 * u[9] - u2 * u[12]);
    du[7] = ζ_hat * (u1 * u[10] - u3 * u[4]);
    du[8] = ζ_hat * (u1 * u[11] - u3 * u[5]);
    du[9] = ζ_hat * (u1 * u[12] - u3 * u[6]);
    du[10] = ζ_hat * (u2 * u[4] - u1 * u[7]);
    du[11] = ζ_hat * (u2 * u[5] - u1 * u[8]);
    du[12] = ζ_hat * (u2 * u[6] - u1 * u[9]);
    # The l1, l2, and l3 do not have a ζ_hat factor because
    # they are already couples per unit initial length;
    # also the case for a distributed body couple due to off-center
    # loading (not considered in this function)
    # [see Moulton et al. Morphoelastic Rods: Part 1]
    du[13] = ζ_hat * (u3 * u[14] - u2 * u[15]) + ζ * n2 + l1; 
    du[14] = ζ_hat * (u1 * u[15] - u3 * u[13]) - ζ * n1 + l2;
    du[15] = ζ_hat * (u2 * u[13] - u1 * u[14]) + l3;
end

## Revise this to consider updated new_K
function self_weight_solve(trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}; 
        m0::Vector{Float64} = [0.0, 0.0, 0.0], uInit::Vector{Float64} = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m0[1], m0[2], m0[3]],
        g::Float64 = -9.8, Ng::Integer = 4, solver = 1, kwargs...) where {T, N}
    activation = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)

    g_range = range(start = 0.0, stop = g, length = Ng)[2:end]
    
    stiffness = trunk.interpolations.K
    ρlin0Int = trunk.ρlin0Int
    u_hat = activation.u_hat

    m10::Float64, m20::Float64, m30::Float64 = m0;
    bcs = SVector{12, Float64}(uInit[1:12]);
    sol = 0;
    Zspan = (0.0, trunk.trunk.L::Float64);

    for gi in g_range
        p = (gi, ρlin0Int, stiffness, u_hat, bcs);
        bvp = TwoPointBVProblem(self_weight_trunk_de!, (self_weight_bc_start!, self_weight_bc_end!), uInit, Zspan, p;
                                bcresid_prototype = (zeros(12), zeros(3)))

        if solver == 1
            # sol = solve(bvp, MIRK4(), kwargs...);
            sol = solve(bvp, MIRK4(), dt = trunk.trunk.L / 100.0, abstol = 1e-3, reltol = 1e-3);
        elseif solver == 2
            # sol = solve(bvp, Shooting(AutoVern7(Rodas4())), dt = trunk.trunk.L / 100.0, abstol = 1e-6, reltol = 1e-6);
        end

        m10, m20, m30 = sol(0)[13:15];
        m_end = sol(trunk.trunk.L)[13:15]
        if (abs(m_end[1]) > 0.001 || abs(m_end[2]) > 0.001 || abs(m_end[3]) > 0.001)
            println("WARNING: Non-zero end moment")
        end
        uInit = [uInit[1], uInit[2], uInit[3], uInit[4], uInit[5], uInit[6], uInit[7], uInit[8], uInit[9], uInit[10], uInit[11], uInit[12], m10, m20, m30];
    end

    sol
end

function build_trunk_ivp(trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}} = ((@SMatrix zeros(3, 5)), (@SMatrix zeros(3, 5)));
        u0::SVector{12, Float64} = SVector{12, Float64}([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
        Zspan = (0.0, trunk.trunk.L),
        kwargs...) where {T, N}
    activation = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
    p = activation.u_hat
    
    ivp = ODEProblem(intrinsic_trunk_de_SA, u0, Zspan, p);
    
    ivp, activation
end

function build_trunk_bvp(trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}} = ((@SMatrix zeros(3, 5)), (@SMatrix zeros(3, 5)));
        bcs = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        m0::Vector{Float64} = [0.0, 0.0, 0.0], uInit::Vector{<:Any} = [bcs..., m0[1], m0[2], m0[3]],
        g::Float64 = -9.8, F::SVector{3, Float64} = SVector{3, Float64}([0.0, 0.0, 0.0])) where {T, N}
    activation = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
    
    stiffness = activation.new_K
    ρlin0Int = trunk.ρlin0Int
    u_hat = activation.u_hat

    bcs = SVector{12, Float64}(bcs)
    Zspan = (0.0, trunk.trunk.L::Float64);
    
    p = (g, ρlin0Int, stiffness, u_hat, F, bcs)
    bvp = TwoPointBVProblem(self_weight_trunk_de!, (self_weight_bc_start!, self_weight_bc_end!), uInit, Zspan, p;
                                bcresid_prototype = (zeros(12), zeros(3)))

    bvp, activation
end

function build_trunk_wrap_nbodyc_bvp(trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}} = ((@SMatrix zeros(3, 5)), (@SMatrix zeros(3, 5)));
        bcs = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        m0::Vector{Float64} = [0.0, 0.0, 0.0], uInit::Vector{Float64} = [bcs..., m0[1], m0[2], m0[3]],
        g::Float64 = -9.8,
        W_C::Float64 = 0.0,
        Z_0::Float64 = trunk.trunk.L) where {T, N}
    activation = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
    L = trunk.trunk.L
    w = W_C / (L - Z_0)

    stiffness = activation.new_K
    ρlin0Int = trunk.ρlin0Int
    u_hat = activation.u_hat

    bcs = SVector{12, Float64}(bcs)
    Zspan = (0.0, L::Float64);

    p = (g, ρlin0Int, stiffness, u_hat, SVector{3, Float64}([w, L, Z_0]), bcs)
    bvp = TwoPointBVProblem(self_weight_wrap_nbodyc_trunk_de!, (self_weight_bc_start!, self_weight_bc_end!), uInit, Zspan, p;
                                bcresid_prototype = (zeros(12), zeros(3)))

    bvp, activation
end


function build_trunk_wrap_bvp(trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}} = ((@SMatrix zeros(3, 5)), (@SMatrix zeros(3, 5)));
        bcs = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        m0::Vector{Float64} = [0.0, 0.0, 0.0], uInit::Vector{Float64} = [bcs..., m0[1], m0[2], m0[3]],
        g::Float64 = -9.8,
        W_C::Float64 = 0.0,
        Z_0::Float64 = trunk.trunk.L,
        r_C::SVector{3, Float64} = SVector{3, Float64}([0.0, 0.0, 0.0])) where {T, N}
    activation = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
    L = trunk.trunk.L
    w = W_C / (L - Z_0)
    R0_act = Z -> (activation.R_factor * trunk.R00 - Z * tan(trunk.trunk.φ0))

    stiffness = activation.new_K
    ρlin0Int = trunk.ρlin0Int
    u_hat = activation.u_hat

    bcs = SVector{12, Float64}(bcs)
    Zspan = (0.0, L::Float64);

    p = (g, ρlin0Int, stiffness, u_hat, (w, L, Z_0, r_C, R0_act), bcs)
    bvp = TwoPointBVProblem(self_weight_wrap_trunk_de!, (self_weight_bc_start!, self_weight_bc_end!), uInit, Zspan, p;
                                bcresid_prototype = (zeros(12), zeros(3)))

    bvp, activation
end

function self_weight_solve_single(bvp::BVProblem, trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}; 
        m0::Vector{Float64} = [0.0, 0.0, 0.0], uInit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m0[1], m0[2], m0[3]], 
        solver = 1, new_bcs = nothing, abstol = 1e-8, reltol = 1e-6, kwargs...) where {T, N}
    a = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
    u_hat = a.u_hat
    
    if isnothing(new_bcs)
        bvp_new = remake(bvp; u0 = uInit, p = (bvp.p[1], bvp.p[2], a.new_K, u_hat, bvp.p[5], bvp.p[6]))
    else
        # println("Remaking BVP bcs")
        bvp_new = remake(bvp; u0 = uInit, p = (bvp.p[1], bvp.p[2], a.new_K, u_hat, bvp.p[5], new_bcs))
        
    end
    # println(bvp_new.p[6])

    if solver == 1
        sol = solve(bvp_new, MIRK4(), dt = trunk.trunk.L / 100.0, abstol = abstol, reltol = reltol);
    elseif solver == 2
        sol = solve(bvp_new, Shooting(AutoVern7(Rodas4())), abstol = abstol, reltol = reltol);
    elseif solver == 3
        sol = solve(bvp_new, BVPM2(), dt = trunk.trunk.L / 100.0, abstol = abstol, reltol = reltol);
    elseif solver == 4
        sol = solve(bvp_new, BVPSOL(), abstol = abstol, reltol = reltol);
    elseif solver == 5
        sol = solve(bvp_new, COLNEW(), dt = trunk.trunk.L / 100.0, abstol = abstol, reltol = reltol);
    end

    m_end = sol(trunk.trunk.L)[13:15]
    # println(m_end)
    # println(sol(0.0)[13:15])
    if (abs(m_end[1]) > 0.001 || abs(m_end[2]) > 0.001 || abs(m_end[3]) > 0.001)
        println("WARNING: Non-zero end moment")
    end

    sol, a
end


# function self_weight_wrap_solve_single(bvp::BVProblem, trunk::TrunkFast{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}; 
#         m0::Vector{Float64} = [0.0, 0.0, 0.0], uInit::Vector{Float64} = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m0[1], m0[2], m0[3]], 
#         solver = 1, new_bcs = nothing, abstol = 1e-4, reltol = 1e-4, kwargs...) where {T, N}
#     a = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)
#     u_hat = a.u_hat
    
#     if isnothing(new_bcs)
#         bvp_new = remake(bvp; u0 = uInit, p = (bvp.p[1], bvp.p[2], a.new_K, u_hat, bvp.p[5], bvp.p[6]))
#     else
#         println("Remaking BVP bcs")
#         bvp_new = remake(bvp; u0 = uInit, p = (bvp.p[1], bvp.p[2], a.new_K, u_hat, bvp.p[5], new_bcs))
#     end
#     println(bvp_new.p[6])

#     if solver == 1
#         # sol = solve(bvp_new, MIRK4(), dt = trunk.trunk.L / 100.0, abstol = 1e-3, reltol = 1e-3);
#         sol = solve(bvp_new, MIRK4(), dt = trunk.trunk.L / 100.0, abstol = abstol, reltol = reltol);
#     elseif solver == 2
#         sol = solve(bvp_new, Shooting(AutoVern7(Rodas4())), dt = trunk.trunk.L / 100.0, abstol = 1e-6, reltol = 1e-6);
#     end

#     m_end = sol(trunk.trunk.L)[13:15]
#     # println(m_end)
#     # println(sol(0.0)[13:15])
#     if (abs(m_end[1]) > 0.001 || abs(m_end[2]) > 0.001 || abs(m_end[3]) > 0.001)
#         println("WARNING: Non-zero end moment")
#     end

#     sol, a
# end


function ivp_solve_single(
        ivp::ODEProblem, 
        trunk::TrunkFast{T, N}, 
        γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}};
        new_u0 = nothing, kwargs...) where {T, N}
    a = ActivatedTrunkQuantities{T, N}(trunkFast = trunk, γ = γ)

    if isnothing(new_u0)
        ivp_new = remake(ivp; p = a.u_hat)
    else
        ivp_new = remake(ivp; u0 = new_u0, p = a.u_hat)
    end

    sol = solve(ivp_new, AutoVern7(Rodas4()), dt = trunk.trunk.L / 100.0, abstol = 1e-8, reltol = 1e-6);

    sol, a
end

function self_weight_bc_start!(residual, u, p)
    bc = p[6];
    residual[1] = u[1] - bc[1];
    residual[2] = u[2] - bc[2];
    residual[3] = u[3] - bc[3];
    residual[4] = u[4] - bc[4];
    residual[5] = u[5] - bc[5];
    residual[6] = u[6] - bc[6];
    residual[7] = u[7] - bc[7];
    residual[8] = u[8] - bc[8];
    residual[9] = u[9] - bc[9];
    residual[10] = u[10] - bc[10];
    residual[11] = u[11] - bc[11];
    residual[12] = u[12] - bc[12];
end

function self_weight_bc_end!(residual, u, p)
    residual[1] = u[13] - 0.0;
    residual[2] = u[14] - 0.0;
    residual[3] = u[15] - 0.0;
end




# function computeH(trunk::Trunk{T, N}, Z_index::Integer, A::Tuple) where {T, N}
#     E = trunk.E
#     δ_lo = deltas_longitudinal(trunk, Z_index, 1)
#     δ_ovo_R = deltas_helical(trunk, Z_index, 2, right)
#     δ_ovo_L = deltas_helical(trunk, Z_index, 2, left)
#     δ_ivo_R = deltas_helical(trunk, Z_index, 3, right)
#     δ_ivo_L = deltas_helical(trunk, Z_index, 3, left)
#     δ_dr = deltas_radial(trunk, Z_index, 4)
#     δ_vr = deltas_radial(trunk, Z_index, 5)
#     AR, AL = A

#     H = E * SVector{3, SVector{N, Float64}}([
#         δ_lo[i] * (AR[i][Z_index, 1] + AL[i][Z_index, 1]) +
#         δ_ovo_R[i] * AR[i][Z_index, 2] +
#         δ_ovo_L[i] * AL[i][Z_index, 2] +  
#         δ_ivo_R[i] * AR[i][Z_index, 3] +
#         δ_ivo_L[i] * AL[i][Z_index, 3] +
#         δ_dr[i] * (AR[i][Z_index, 4] + AL[i][Z_index, 4]) +
#         δ_vr[i] * (AR[i][Z_index, 5] + AL[i][Z_index, 5])
#         for i in 1:3
#             ])
    
#     H
# end

# function computeH(trunk::Trunk{T, N}, γ::Tuple{SMatrix{T, 5, Float64}, SMatrix{T, 5, Float64}}) where {T, N}
#     A = constant_activation(trunk, γ)
    
#     H = SVector{3, SMatrix{T, N, Float64}}
#     H = SVector{T, SVector{3, SVector{N, Float64}}}
#     for i in 1:T
#         H = computeH(trunk, i, A) # this is an SVector{3, SVector{N, Float64}}
        
        
#         println(typeof(H))
#         println(H)
#     end


# end