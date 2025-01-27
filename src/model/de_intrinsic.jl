######### Kirchhoff rod kinematics and mechanics
######### Without external forces
#region ===========================
"""
    $(TYPEDSIGNATURES)

In-place DE function for the kinematics of the Kirchhoff filament
for symbolic inputs.

The parameter `p` stores the runtime-generated functions for
`ζ_hat`, `u1_hat`, `u2_hat`, and `u3_hat`.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`
"""
function intrinsicConfDESym!(du, u, p, Z)
    ζ_hat, u1_hat, u2_hat, u3_hat = [p[1](Z), p[2](Z), p[3](Z), p[4](Z)]
    du[1] = ζ_hat * u[10]
    du[2] = ζ_hat * u[11]
    du[3] = ζ_hat * u[12]
    du[4] = ζ_hat * (u3_hat * u[7] - u2_hat * u[10])
    du[5] = ζ_hat * (u3_hat * u[8] - u2_hat * u[11])
    du[6] = ζ_hat * (u3_hat * u[9] - u2_hat * u[12])
    du[7] = ζ_hat * (u1_hat * u[10] - u3_hat * u[4])
    du[8] = ζ_hat * (u1_hat * u[11] - u3_hat * u[5])
    du[9] = ζ_hat * (u1_hat * u[12] - u3_hat * u[6])
    du[10] = ζ_hat * (u2_hat * u[4] - u1_hat * u[7])
    du[11] = ζ_hat * (u2_hat * u[5] - u1_hat * u[8])
    du[12] = ζ_hat * (u2_hat * u[6] - u1_hat * u[9])
end

"""
    $(TYPEDSIGNATURES)

DE function for the kinematics of the Kirchhoff filament
for symbolic inputs.

The parameter `p` stores the runtime-generated functions for
`ζ_hat`, `u1_hat`, `u2_hat`, and `u3_hat`.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

Note: this function is NOT in-place.
"""
function intrinsicConfDESym(u, p, Z)
    ζ_hat, u1_hat, u2_hat, u3_hat = [p[1](Z), p[2](Z), p[3](Z), p[4](Z)]
    du1 = ζ_hat * u[10]
    du2 = ζ_hat * u[11]
    du3 = ζ_hat * u[12]
    du4 = ζ_hat * (u3_hat * u[7] - u2_hat * u[10])
    du5 = ζ_hat * (u3_hat * u[8] - u2_hat * u[11])
    du6 = ζ_hat * (u3_hat * u[9] - u2_hat * u[12])
    du7 = ζ_hat * (u1_hat * u[10] - u3_hat * u[4])
    du8 = ζ_hat * (u1_hat * u[11] - u3_hat * u[5])
    du9 = ζ_hat * (u1_hat * u[12] - u3_hat * u[6])
    du10 = ζ_hat * (u2_hat * u[4] - u1_hat * u[7])
    du11 = ζ_hat * (u2_hat * u[5] - u1_hat * u[8])
    du12 = ζ_hat * (u2_hat * u[6] - u1_hat * u[9])
    SVector{12}(du1, du2, du3, du4, du5, du6, du7, du8, du9, du10, du11, du12)
end

"""
    $(TYPEDSIGNATURES)

DE function for the kinematics of the Kirchhoff filament
for numerical Float32 inputs and outputs. To be used with GPU computations.

Use a `StaticArray` for the `u` input.

The parameter `p` stores precomputed quantities in an `SMatrix`.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

The output is an `SVector` array to accelerate computations.

Note: this function is NOT in-place.
"""
function intrinsicConfDEGPU(u, p, Z)
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHatGPU(Z, p)
    du1 = ζ_hat * u[10]
    du2 = ζ_hat * u[11]
    du3 = ζ_hat * u[12]
    du4 = ζ_hat * (u3_hat * u[7] - u2_hat * u[10])
    du5 = ζ_hat * (u3_hat * u[8] - u2_hat * u[11])
    du6 = ζ_hat * (u3_hat * u[9] - u2_hat * u[12])
    du7 = ζ_hat * (u1_hat * u[10] - u3_hat * u[4])
    du8 = ζ_hat * (u1_hat * u[11] - u3_hat * u[5])
    du9 = ζ_hat * (u1_hat * u[12] - u3_hat * u[6])
    du10 = ζ_hat * (u2_hat * u[4] - u1_hat * u[7])
    du11 = ζ_hat * (u2_hat * u[5] - u1_hat * u[8])
    du12 = ζ_hat * (u2_hat * u[6] - u1_hat * u[9])
    SVector{12}(du1, du2, du3, du4, du5, du6, du7, du8, du9, du10, du11, du12)
end

"""
    $(TYPEDSIGNATURES)

DE function for the kinematics of the Kirchhoff filament
for numerical inputs. The function uses `StaticArray`
wherever possible.

Use a `StaticArray` for the `u` input.

The parameter `p` stores precomputed quantities in an `SMatrix`.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

The output is an `SVector` array to accelerate computations.

Note: this function is NOT in-place.
"""
function intrinsicConfDESA(u, p, Z)
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHat(Z, p)
    du1 = ζ_hat * u[10]
    du2 = ζ_hat * u[11]
    du3 = ζ_hat * u[12]
    du4 = ζ_hat * (u3_hat * u[7] - u2_hat * u[10])
    du5 = ζ_hat * (u3_hat * u[8] - u2_hat * u[11])
    du6 = ζ_hat * (u3_hat * u[9] - u2_hat * u[12])
    du7 = ζ_hat * (u1_hat * u[10] - u3_hat * u[4])
    du8 = ζ_hat * (u1_hat * u[11] - u3_hat * u[5])
    du9 = ζ_hat * (u1_hat * u[12] - u3_hat * u[6])
    du10 = ζ_hat * (u2_hat * u[4] - u1_hat * u[7])
    du11 = ζ_hat * (u2_hat * u[5] - u1_hat * u[8])
    du12 = ζ_hat * (u2_hat * u[6] - u1_hat * u[9])
    SVector{12}(du1, du2, du3, du4, du5, du6, du7, du8, du9, du10, du11, du12)
end

"""
    $(TYPEDSIGNATURES)

DE function for the kinematics of the Kirchhoff filament
for numerical inputs.

The parameter `p` stores precomputed quantities
in `PrecomputedQuantities`.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`
"""
function intrinsicConfDE!(du, u, p, Z)
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHatSym(Z, p)
    du[1:3] = ζ_hat * u[10:12]
    du[4:6] = ζ_hat * (u3_hat * u[7:9] - u2_hat * u[10:12])
    du[7:9] = ζ_hat * (u1_hat * u[10:12] - u3_hat * u[4:6])
    du[10:12] = ζ_hat * (u2_hat * u[4:6] - u1_hat * u[7:9])
end

"""
    $(TYPEDSIGNATURES)

In-place DE function for the kinematics of the Kirchhoff filament
for numerical inputs.

The parameter `p` stores precomputed quantities
in `PrecomputedQuantities`.

The unknown function u has the following structure:
`(u[1], u[2], u[3])` = r = filament centerline along `Z`

`(
    (u[4], u[5], u[6]),
    (u[7], u[8], u[9]),
    (u[10], u[11], u[12])
)`   = `(d1, d2, d3)` = director basis along `Z`

Note: this function is NOT in-place.
"""
function intrinsicConfDE(u, p, Z)
    ζ_hat, u1_hat, u2_hat, u3_hat = computeUHatSym(Z, p)
    du = MVector{12, Float64}(zeros(12))
    du[1:3] = ζ_hat * u[10:12]
    du[4:6] = ζ_hat * (u3_hat * u[7:9] - u2_hat * u[10:12])
    du[7:9] = ζ_hat * (u1_hat * u[10:12] - u3_hat * u[4:6])
    du[10:12] = ζ_hat * (u2_hat * u[4:6] - u1_hat * u[7:9])
    SVector{12}(du)
end
#endregion ===========================
