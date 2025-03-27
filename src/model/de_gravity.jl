############################################################################
### Kirchhoff rod kinematics and mechanics:                              ###
### Computation of deformed filament configurations with external forces ###
############################################################################
#region ===========================
"""
    $(TYPEDSIGNATURES)

In-place DE function for the mechanics of the Kirchhoff filament
under gravity for non-symbolic computations.

The parameter `p` has the following structure:
-   `p[1]` = `g::Float64` = gravitational acceleration
-   `p[2]` = `ρlinInt::Function` = linear density integral function
-   `p[3]` = `SVector{4, Float64}` = stiffness vector
-   `p[4]` = `precomputedQuantities::SMatrix` = precomputed quantities matrix

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

`(u[13], u[14], u[15])` = moments along `Z`
"""
function selfWeightDE!(
        du,
        u,
        p::Tuple{
            Float64, Function, SVector{4, Float64}, <:AbstractMatrix, SVector{12, Float64}},
        Z
)
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHat(Z, p[4])

    ρlinInt = p[2](Z)

    n1 = -p[1] * ρlinInt * u[6]
    n2 = -p[1] * ρlinInt * u[9]
    n3 = -p[1] * ρlinInt * u[12]

    ζF = n3 / p[3][1] + 1.0
    ζ = ζ_hat * ζF
    u1 = u1_hat + u[13] / p[3][2]
    u2 = u2_hat + u[14] / p[3][3]
    u3 = u3_hat + u[15] / p[3][4]

    du[1] = ζ * u[10]
    du[2] = ζ * u[11]
    du[3] = ζ * u[12]
    du[4] = ζ_hat * (u3 * u[7] - u2 * u[10])
    du[5] = ζ_hat * (u3 * u[8] - u2 * u[11])
    du[6] = ζ_hat * (u3 * u[9] - u2 * u[12])
    du[7] = ζ_hat * (u1 * u[10] - u3 * u[4])
    du[8] = ζ_hat * (u1 * u[11] - u3 * u[5])
    du[9] = ζ_hat * (u1 * u[12] - u3 * u[6])
    du[10] = ζ_hat * (u2 * u[4] - u1 * u[7])
    du[11] = ζ_hat * (u2 * u[5] - u1 * u[8])
    du[12] = ζ_hat * (u2 * u[6] - u1 * u[9])
    du[13] = ζ_hat * (u3 * u[14] - u2 * u[15]) + ζ * n2
    du[14] = ζ_hat * (u1 * u[15] - u3 * u[13]) - ζ * n1
    du[15] = ζ_hat * (u2 * u[13] - u1 * u[14])
end

# Tapered
function selfWeightDE!(
        du,
        u,
        p::Tuple{Float64, Function, SVector{4, <:AbstractInterpolation}, Tuple},
        Z
)
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHat(Z, p[4])
    ρlinInt = p[2](Z)

    n1 = -p[1] * ρlinInt * u[6]
    n2 = -p[1] * ρlinInt * u[9]
    n3 = -p[1] * ρlinInt * u[12]

    ζF = n3 / p[3][1](Z) + 1.0
    ζ = ζ_hat * ζF
    u1 = u1_hat + u[13] / p[3][2](Z)
    u2 = u2_hat + u[14] / p[3][3](Z)
    u3 = u3_hat + u[15] / p[3][4](Z)

    du[1] = ζ * u[10]
    du[2] = ζ * u[11]
    du[3] = ζ * u[12]
    du[4] = ζ_hat * (u3 * u[7] - u2 * u[10])
    du[5] = ζ_hat * (u3 * u[8] - u2 * u[11])
    du[6] = ζ_hat * (u3 * u[9] - u2 * u[12])
    du[7] = ζ_hat * (u1 * u[10] - u3 * u[4])
    du[8] = ζ_hat * (u1 * u[11] - u3 * u[5])
    du[9] = ζ_hat * (u1 * u[12] - u3 * u[6])
    du[10] = ζ_hat * (u2 * u[4] - u1 * u[7])
    du[11] = ζ_hat * (u2 * u[5] - u1 * u[8])
    du[12] = ζ_hat * (u2 * u[6] - u1 * u[9])
    du[13] = ζ_hat * (u3 * u[14] - u2 * u[15]) + ζ * n2
    du[14] = ζ_hat * (u1 * u[15] - u3 * u[13]) - ζ * n1
    du[15] = ζ_hat * (u2 * u[13] - u1 * u[14])
end

"""
    $(TYPEDSIGNATURES)

In-place DE function for the mechanics of the Kirchhoff filament
under gravity for symbolic computations.

The parameter `p` has the following structure:
-   `p[1]` = `g::Float64` = gravitational acceleration
-   `p[2]` = `filament::AFilament`
-   `p[3]` = `u_f` = runtime-generated functions for `ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`, and `ζ_hat_AD`

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

`(u[13], u[14], u[15])` = moments along `Z`
"""
function selfWeightDESym!(du, u, p, Z)
    ### p[1] = g
    ### p[2] = filament
    ### p[3] = u_f
    u_f = p[3]
    ζ_hat, u1_hat, u2_hat, u3_hat = [u_f[1](Z), u_f[2](Z), u_f[3](Z), u_f[4](Z)]
    filament = p[2]

    # NOTE 1: For constant zetaHat, this is the same as the old version ρlinInt = (filament.m / p[4]) * evaluate_integral_AD(u_f[5], Z, filament.L)
    # NOTE 2: For piecewise constant helical angle (hence, piecewise constant zetaHat), the old version
    # ρlinInt = (filament.m / p[4]) * evaluate_integral_AD(u_f[5], Z, filament.L)
    # is quite close numerically (difference impercetible in plots) to the corrected version below. This happens if the zetaHat in the 
    # helical portion is not significantly different than the zetaHat in the longitudinal portion or if the helical portion is a 
    # small segment of the overall filament.
    ρlinInt = filament.m / filament.L * (filament.L - Z) # Assuming homogeneous volumetric density in initial conf. and a non-tapered filament

    n1 = -p[1] * ρlinInt * u[6]
    n2 = -p[1] * ρlinInt * u[9]
    n3 = -p[1] * ρlinInt * u[12]

    ζF = n3 / filament.stiffness.K0 + 1.0
    ζ = ζ_hat * ζF
    u1 = u1_hat + u[13] / filament.stiffness.K1
    u2 = u2_hat + u[14] / filament.stiffness.K2
    u3 = u3_hat + u[15] / filament.stiffness.K3

    du[1] = ζ * u[10]
    du[2] = ζ * u[11]
    du[3] = ζ * u[12]
    du[4] = ζ_hat * (u3 * u[7] - u2 * u[10])
    du[5] = ζ_hat * (u3 * u[8] - u2 * u[11])
    du[6] = ζ_hat * (u3 * u[9] - u2 * u[12])
    du[7] = ζ_hat * (u1 * u[10] - u3 * u[4])
    du[8] = ζ_hat * (u1 * u[11] - u3 * u[5])
    du[9] = ζ_hat * (u1 * u[12] - u3 * u[6])
    du[10] = ζ_hat * (u2 * u[4] - u1 * u[7])
    du[11] = ζ_hat * (u2 * u[5] - u1 * u[8])
    du[12] = ζ_hat * (u2 * u[6] - u1 * u[9])
    du[13] = ζ_hat * (u3 * u[14] - u2 * u[15]) + ζ * n2
    du[14] = ζ_hat * (u1 * u[15] - u3 * u[13]) - ζ * n1
    du[15] = ζ_hat * (u2 * u[13] - u1 * u[14])
end

"""
    $(TYPEDSIGNATURES)

In-place DE function for the mechanics of the Kirchhoff filament
under gravity. Use for non-symbolic computation on the GPU.

The parameter `p` has the following structure:
-   `p[1]` = `g::Float32` = gravitational acceleration
-   `p[2]` = `filament::AFilamentGPU`
-   `p[3]` = `precomputedQuantities::SMatrix` = precomputed quantities matrix

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

`(u[13], u[14], u[15])` = moments along `Z`
"""
function selfWeightDE_GPU!(du, u, p, Z) # has to be in-place for BVP solvers?
    # p[1] = g
    # p[2] = filament
    # p[3] = precomputedQuantities
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHatGPU(Z, p[3])
    filament = p[2]

    ρlinInt = ζ_hat * filament.ρlin * (filament.L - Z) # Assuming constant ζ_hat

    n1 = -p[1] * ρlinInt * u[6]
    n2 = -p[1] * ρlinInt * u[9]
    n3 = -p[1] * ρlinInt * u[12]

    ζF = n3 / filament.stiffness.K0 + 1.0f0
    ζ = ζ_hat * ζF
    u1 = u1_hat + u[13] / filament.stiffness.K1
    u2 = u2_hat + u[14] / filament.stiffness.K2
    u3 = u3_hat + u[15] / filament.stiffness.K3

    du[1] = ζ * u[10]
    du[2] = ζ * u[11]
    du[3] = ζ * u[12]
    du[4] = ζ * (u3 * u[7] - u2 * u[10])
    du[5] = ζ * (u3 * u[8] - u2 * u[11])
    du[6] = ζ * (u3 * u[9] - u2 * u[12])
    du[7] = ζ * (u1 * u[10] - u3 * u[4])
    du[8] = ζ * (u1 * u[11] - u3 * u[5])
    du[9] = ζ * (u1 * u[12] - u3 * u[6])
    du[10] = ζ * (u2 * u[4] - u1 * u[7])
    du[11] = ζ * (u2 * u[5] - u1 * u[8])
    du[12] = ζ * (u2 * u[6] - u1 * u[9])
    du[13] = u3 * u[14] - u2 * u[15] + n2
    du[14] = u1 * u[15] - u3 * u[13] - n1
    du[15] = u2 * u[13] - u1 * u[14]
end

"""
    $(TYPEDSIGNATURES)

In-place boundary condition function that defines the
residuals that are to become identically zero to solve the
corresponding BVP. To be used with numerical Float32 GPU computation.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

`(u[13], u[14], u[15])` = moments along `Z`

The form of the residuals implies that:
-   `r(0)` = (0, 0, 0)
-   `d1(0)` = (1, 0, 0)
-   `d2(0)` = (0, 1, 0)
-   `d3(0)` = (0, 0, 1)
-   `m(L)`  = (0, 0, 0)
"""
function selfWeightBC_GPU!(residual, u, p, t)
    residual[1] = u[1][1] - 0.0f0
    residual[2] = u[1][2] - 0.0f0
    residual[3] = u[1][3] - 0.0f0
    residual[4] = u[1][4] - 1.0f0
    residual[5] = u[1][5] - 0.0f0
    residual[6] = u[1][6] - 0.0f0
    residual[7] = u[1][7] - 0.0f0
    residual[8] = u[1][8] - 1.0f0
    residual[9] = u[1][9] - 0.0f0
    residual[10] = u[1][10] - 0.0f0
    residual[11] = u[1][11] - 0.0f0
    residual[12] = u[1][12] - 1.0f0
    residual[13] = u[end][13] - 0.0f0
    residual[14] = u[end][14] - 0.0f0
    residual[15] = u[end][15] - 0.0f0
end

"""
    $(TYPEDSIGNATURES)

In-place boundary condition function that defines the
residuals that are to become identically zero to solve the
corresponding BVP. To be used with symbolic or numerical
computations.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

`(u[13], u[14], u[15])` = moments along `Z`

The form of the residuals implies that:
-   `r(0)` = p[5][1:3]
-   `d1(0)` = p[5][4:6]
-   `d2(0)` = p[5][7:9]
-   `d3(0)` = p[5][10:12]
-   `m(L)`  = (0, 0, 0)
"""
function selfWeightBC!(residual, u, p, t)
    bc = p[5]
    residual[1] = u[1][1] - bc[1]
    residual[2] = u[1][2] - bc[2]
    residual[3] = u[1][3] - bc[3]
    residual[4] = u[1][4] - bc[4]
    residual[5] = u[1][5] - bc[5]
    residual[6] = u[1][6] - bc[6]
    residual[7] = u[1][7] - bc[7]
    residual[8] = u[1][8] - bc[8]
    residual[9] = u[1][9] - bc[9]
    residual[10] = u[1][10] - bc[10]
    residual[11] = u[1][11] - bc[11]
    residual[12] = u[1][12] - bc[12]
    residual[13] = u[end][13] - 0.0
    residual[14] = u[end][14] - 0.0
    residual[15] = u[end][15] - 0.0
end

"""
    $(TYPEDSIGNATURES)

In-place boundary condition function that defines the
residuals that are to become identically zero at Z = 0 to solve the
corresponding two-point BVP.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

The form of the residuals implies that:
-   `r(0)` = p[5][1:3]
-   `d1(0)` = p[5][4:6]
-   `d2(0)` = p[5][7:9]
-   `d3(0)` = p[5][10:12]
"""
function selfWeightBCStart!(residual, u, p)
    bc = p[5]
    residual[1] = u[1] - bc[1]
    residual[2] = u[2] - bc[2]
    residual[3] = u[3] - bc[3]
    residual[4] = u[4] - bc[4]
    residual[5] = u[5] - bc[5]
    residual[6] = u[6] - bc[6]
    residual[7] = u[7] - bc[7]
    residual[8] = u[8] - bc[8]
    residual[9] = u[9] - bc[9]
    residual[10] = u[10] - bc[10]
    residual[11] = u[11] - bc[11]
    residual[12] = u[12] - bc[12]
end

"""
    $(TYPEDSIGNATURES)

In-place boundary condition function that defines the
residuals that are to become identically zero at Z = L to solve the
corresponding two-point BVP.

The unknown function u has the following structure:
`(u[13], u[14], u[15])` = moments along `Z`
"""
function selfWeightBCEnd!(residual, u, p)
    residual[1] = u[13] - 0.0
    residual[2] = u[14] - 0.0
    residual[3] = u[15] - 0.0
end
#endregion ===========================
