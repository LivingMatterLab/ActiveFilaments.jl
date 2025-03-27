#####################################################################
### Mathematical functions relevant to the Active Filament theory ###
#####################################################################
#region ===========================
"""
    $(TYPEDSIGNATURES)

Computes the stiffnesses `K0`, `K1`, `K2`, and `K3` for a given filament `ring`.
"""
function computeK(ring::Ring)
    ν = ring.mechanicalProperties.ν
    E = ring.mechanicalProperties.E
    R1 = ring.geometry.R1
    R2 = ring.geometry.R2

    K0 = E * pi * (R2 .^ 2 - R1 .^ 2)
    K1 = E * pi / 4 * (R2 .^ 4 - R1 .^ 4)
    K2 = K1
    K3 = E * pi / (4 * (1 + ν)) * (R2 .^ 4 - R1 .^ 4)

    SVector{4}(K0, K1, K2, K3)
end

"""
    $(TYPEDSIGNATURES)

Computes the stiffnesses `K0`, `K1`, `K2`, and `K3` for the inner tube
`InnerTube` of the filament.
"""
function computeK(innerTube::InnerTube)
    ν = innerTube.mechanicalProperties.ν
    E = innerTube.mechanicalProperties.E
    R1 = innerTube.geometry.R1
    R2 = innerTube.geometry.R2

    K0 = E * pi * (R2 .^ 2 - R1 .^ 2)
    K1 = E * pi / 4 * (R2 .^ 4 - R1 .^ 4)
    K2 = K1
    K3 = E * pi / (4 * (1 + ν)) * (R2 .^ 4 - R1 .^ 4)

    SVector{4}(K0, K1, K2, K3)
end

"""
    $(TYPEDSIGNATURES)

Computes the stiffnesses `K0`, `K1`, `K2`, and `K3`
for the entire filament.
"""
function computeK(rings::Vector{Ring{T}} where {T}, innerTube::InnerTube{T} where {T})
    K = computeK(innerTube)
    for ring in rings
        K_ring = computeK(ring)
        K += K_ring
    end
    K
end

"""
    $(TYPEDSIGNATURES)

Computes the delta quantities `δ0`, `δ1`, `δ2`, `δ3` for the filament `ring`.
"""
function computeDeltas(ring::Ring)
    α2_st = ring.fiberArchitecture.α2
    R1 = ring.geometry.R1
    R2 = ring.geometry.R2

    if α2_st isa PiecewiseStructure # Check if helical angle is piecewise in Z
        expressionsδ0 = []
        expressionsδ1 = []
        expressionsδ2 = []
        expressionsδ3 = []
        for (α2, range) in zip(α2_st.expressions, α2_st.piecewiseRanges)
            if (α2 !== 0 && α2 !== 0.0f0 && α2 !== 0.0) # Helical
                ν = ring.mechanicalProperties.ν
                δ0 = (ν + 1) * R2^2 * cot(α2)^2 *
                     log(R2^2 * sec(α2)^2 / (R1^2 * tan(α2)^2 + R2^2)) + ν * (R1^2 - R2^2)
                δ1 = (R1 - R2) *
                     (ν * (R1^2 + R2 * R1 + R2^2) - 3 * (ν + 1) * R2^2 * cot(α2)^2) +
                     3 * (ν + 1) * R2^3 * cot(α2)^3 * atan(R1 * tan(α2) / R2) -
                     3 * α2 * (ν + 1) * R2^3 * cot(α2)^3
                δ2 = δ1
                δ3 = -(1 + ν) * R2 * cot(α2) *
                     (
                         R1^2 - R2^2 +
                         R2^2 * cot(α2)^2 *
                         log(2 * R2^2 / ((R2^2 - R1^2) * cos(2 * α2) + R1^2 + R2^2))
                     )
                push!(expressionsδ0, δ0)
                push!(expressionsδ1, δ1)
                push!(expressionsδ2, δ2)
                push!(expressionsδ3, δ3)
            else # Longitudinal
                δ0 = R2^2 - R1^2
                δ1 = R2^3 - R1^3
                δ2 = δ1
                δ3 = 0
                push!(expressionsδ0, δ0)
                push!(expressionsδ1, δ1)
                push!(expressionsδ2, δ2)
                push!(expressionsδ3, δ3)
            end
        end
        δ0 = PiecewiseStructure(expressionsδ0, α2_st.piecewiseRanges)
        δ1 = PiecewiseStructure(expressionsδ1, α2_st.piecewiseRanges)
        δ2 = PiecewiseStructure(expressionsδ2, α2_st.piecewiseRanges)
        δ3 = PiecewiseStructure(expressionsδ3, α2_st.piecewiseRanges)
    else # Helical angle constant in Z
        if ring.geometry.phi2 == 0 # Non-tapered
            α2 = α2_st
            if (α2 !== 0 && α2 !== 0.0f0 && α2 !== 0.0 && α2 !== -0.0) # Helical case
                ν = ring.mechanicalProperties.ν
                δ0 = (ν + 1) * R2^2 * cot(α2)^2 *
                     log(R2^2 * sec(α2)^2 / (R1^2 * tan(α2)^2 + R2^2)) + ν * (R1^2 - R2^2)
                δ1 = (R1 - R2) *
                     (ν * (R1^2 + R2 * R1 + R2^2) - 3 * (ν + 1) * R2^2 * cot(α2)^2) +
                     3 * (ν + 1) * R2^3 * cot(α2)^3 * atan(R1 * tan(α2) / R2) -
                     3 * α2 * (ν + 1) * R2^3 * cot(α2)^3
                δ2 = δ1
                δ3 = -(1 + ν) * R2 * cot(α2) *
                     (
                         R1^2 - R2^2 +
                         R2^2 * cot(α2)^2 *
                         log(2 * R2^2 / ((R2^2 - R1^2) * cos(2 * α2) + R1^2 + R2^2))
                     )
            else # Longitudinal case
                δ0 = R2^2 - R1^2
                δ1 = R2^3 - R1^3
                δ2 = δ1
                δ3 = 0
            end
        else # Tapered
            α2 = α2_st
            cph = tan(ring.geometry.phi2) ./ R2
            ν = ring.mechanicalProperties.ν

            cph2 = cph .^ 2
            R12 = R1 .^ 2
            R22 = R2 .^ 2
            R13 = R12 .* R1
            R23 = R22 .* R2
            R1cph = R1 .* cph
            R2cph = R2 .* cph
            tph1 = 1 .+ R12 .* cph2
            tph2 = 1 .+ R22 .* cph2

            if (α2 !== 0 && α2 !== 0.0f0 && α2 !== 0.0 && α2 !== -0.0) # Helical
                ca = tan(α2) ./ R2
                ca2 = ca .^ 2
                R1ca = R1 .* ca
                R2ca = R2 .* ca
                ta1 = 1 .+ R12 .* ca2
                ta2 = 1 .+ R22 .* ca2
                cpha2 = cph2 - ca2
                sq1 = sqrt.(tph1 .* Complex.(cpha2))
                sq2 = sqrt.(tph2 .* Complex.(cpha2))
                R1cph2 = R1 .* cph2
                R2cph2 = R2 .* cph2

                δ0 = 2 * (R12 - R22) * ν -
                     (2 * (1 + ν) * log.((tph1 .* ta2) ./ (tph2 .* ta1))) ./ (cpha2)

                δ1 = (
                    2 ./ (cph .* ca .* (cpha2))
                    .*
                    (
                    3 * (1 + ν) * (atan.(R1cph) - atan.(R2cph)) .* ca
                    + (R13 - R23) * ν .* (cph2 .* cph) .* ca
                    +
                    cph .* (
                        -3 * (1 + ν) * (atan.(R1ca) - atan.(R2ca)) +
                        (R23 - R13) * ν .* (ca2 .* ca)
                    )
                )
                )
                δ2 = δ1
                δ3 = (
                    1 ./ (cph2 .* ca2 .* sqrt.(Complex.(cpha2)))
                    .*
                    (
                    3 * cph2 .* (
                    atan.((im * R1cph2 + ca) ./ sq1)
                    -
                    atan.((im * R2cph2 + ca) ./ sq2)
                    -
                    im * atanh.((R1cph2 + im * ca) ./ sq1)
                    +
                    im * atanh.((R2cph2 + im * ca) ./ sq2)
                )
                    +
                    6 * ca .* (sq2 - sq1)
                )
                )
                return SVector{4}(real(δ0), real(δ1), real(δ2), real(δ3))
            else # Longitudinal
                δ0 = 2 * (R12 - R22) * ν - (2 * (1 + ν) * log.(tph1 ./ tph2)) ./ cph2
                δ1 = 2 * (R13 - R23) * ν +
                     (6 * (1 + ν) * (atan.(R1cph) - atan.(R2cph) + (R2 - R1) .* cph)) ./
                     (cph2 .* cph)
                δ2 = δ1
                δ3 = 0
            end
        end
    end

    return SVector{4}(δ0, δ1, δ2, δ3)
end

"""
    $(TYPEDSIGNATURES)

Computes the symbolic or numerical prefactors 
for the `H0`, `H1`, `H2`, and `H3` expressions in a given filament `ring`.
"""
function computeHPrefactors(ring::Ring)
    ν = ring.mechanicalProperties.ν
    E = ring.mechanicalProperties.E

    δ = computeDeltas(ring)

    if δ[1] isa PiecewiseStructure # Check if deltas are piecewise-constant in Z
        pre_H0 = simplify(E * pi / 2 * δ[1])
        pre_H1 = simplify(-E * pi / 3 * δ[2])
        pre_H2 = simplify(E * pi / 3 * δ[3])
        pre_H3 = simplify(E * pi / 2 * δ[4] / (1 + ν))
    else # Deltas not piecewise-constant
        if (ring.geometry.phi2 == 0) # Non-tapered
            pre_H0 = E * pi / 2 * δ[1]
            pre_H1 = -E * pi / 3 * δ[2]
            pre_H2 = E * pi / 3 * δ[3]
            pre_H3 = E * pi / 2 * δ[4] / (1 + ν)
        else # Tapered
            pre_H0 = E * pi / 4 * δ[1]
            pre_H1 = -E * pi / 6 * δ[2]
            pre_H2 = E * pi / 6 * δ[3]
            pre_H3 = E * pi / 6 * δ[4]
        end
    end

    return SVector{4}(pre_H0, pre_H1, pre_H2, pre_H3)
end

"""
    $(TYPEDSIGNATURES)

Computes the numerical prefactors 
for the (`ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`) quantities.  
"""
function computePropertyPrefactors(filament::AFilament)
    stiffness = filament.stiffness
    K0 = stiffness.K0
    K1 = stiffness.K1
    K2 = stiffness.K2
    K3 = stiffness.K3

    prefactors = Vector{Prefactors}()
    for ring in filament.rings
        pre_H = computeHPrefactors(ring)
        push!(
            prefactors,
            Prefactors(
                simplify(pre_H[1] ./ K0),
                simplify(pre_H[2] ./ K1),
                simplify(pre_H[3] ./ K2),
                simplify(pre_H[4] ./ K3)
            )
        )
    end

    prefactors
end

"""
    $(TYPEDSIGNATURES)

Computes the symbolic prefactors 
for the (`ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`) quantities.  
"""
function computePropertyPrefactorsSym(filament::AFilament)
    stiffness = filament.stiffness
    K0 = stiffness.K0
    K1 = stiffness.K1
    K2 = stiffness.K2
    K3 = stiffness.K3

    prefactors = Vector{Prefactors}()
    for ring in filament.rings
        pre_H = computeHPrefactors(ring)
        push!(
            prefactors,
            Prefactors(
                simplify(pre_H[1] / K0),
                simplify(pre_H[2] / K1),
                simplify(pre_H[3] / K2),
                simplify(pre_H[4] / K3)
            )
        )
    end

    prefactors
end

"""
    $(TYPEDSIGNATURES)

Converts a piecewise angular activation defined by the `ActivationPiecewiseGamma`
input to the equivalent `ActivationFourier`.
"""
function piecewiseGammaToFourier(activationGamma::ActivationPiecewiseGamma)
    γ = activationGamma.γ
    σ = activationGamma.σ
    θ0 = activationGamma.θ0
    N = length(γ)

    a0 = σ / pi * sum(γ)
    a1 = 0.0
    b1 = 0.0
    for i in 0:(N - 1)
        a1 += γ[i + 1] * cos(θ0 + 2 * pi * i / N)
        b1 += γ[i + 1] * sin(θ0 + 2 * pi * i / N)
    end
    a1 *= 2 * sin(σ / 2) / pi
    b1 *= 2 * sin(σ / 2) / pi

    ActivationFourier(a0, a1, b1)
end

"""
    $(TYPEDSIGNATURES)

Computes the `PrecomputedQuantities` for a given set of `ActivationFourier`
activations in the `filament`.

The resulting `PrecomputedQuantities` are used to speed up computation.
"""
function computeUQuantities(
        filament::AFilament{1, M} where {M},
        activationsFourier::Vector{ActivationFourier},
        prefactors::Vector{Prefactors};
        interp = cubic_spline_interpolation
)
    Z = filament.Z
    M = typeof(filament).parameters[2]

    uPrefactors = Vector{Prefactors{Interpolations.Extrapolation}}()
    ϕ = Vector{Float64}()
    argTerms = Vector{Interpolations.Extrapolation}()

    for (activationFourier, propertyPrefactors, ring) in zip(
        activationsFourier, prefactors, filament.rings)
        a0 = activationFourier.a0
        a1 = activationFourier.a1
        b1 = activationFourier.b1
        A = sqrt(a1^2 + b1^2)
        R2 = ring.geometry.R2
        α2 = ring.fiberArchitecture.α2
        push!(ϕ, -atan(b1, a1))

        # Theta2 solution for the linear tapering case
        theta2T = -tan(α2) * log.(R2 / R2[1]) / sin(ring.geometry.phi2)
        theta2T_interp = interp(filament.Z, theta2T)
        push!(argTerms, theta2T_interp)

        push!(
            uPrefactors,
            Prefactors{Interpolations.Extrapolation}(
                interp(Z, propertyPrefactors.pre_ζ * a0),
                interp(Z, propertyPrefactors.pre_u1 * A),
                interp(Z, -propertyPrefactors.pre_u2 * A),
                interp(Z, propertyPrefactors.pre_u3 * a0)
            )
        )
    end

    PrecomputedQuantities{Float64, Interpolations.Extrapolation}(uPrefactors, ϕ, argTerms)
end

### The functions computeUQuantities and computeUQuantitiesSym can be combined 
### into one with other necessary changes, because they are the same apart from 
### the ActivationFourier type in the collection
function computeUQuantities(
        filament::AFilament{0, M} where {M},
        activationsFourier::Vector{ActivationFourier},
        prefactors::Vector{Prefactors}
)
    uPrefactors = Vector{Prefactors{Float64}}()
    ϕ = Vector{Float64}()
    argTerms = Vector{Float64}()
    for (activationFourier, propertyPrefactors, ring) in zip(
        activationsFourier, prefactors, filament.rings)
        a0 = activationFourier.a0
        a1 = activationFourier.a1
        b1 = activationFourier.b1
        A = sqrt(a1^2 + b1^2)
        R2 = ring.geometry.R2
        α2 = ring.fiberArchitecture.α2
        push!(ϕ, simplify(-atan(b1, a1)))
        push!(argTerms, simplify(tan(α2) / R2))

        push!(
            uPrefactors,
            Prefactors{Float64}(
                simplify(propertyPrefactors.pre_ζ * a0),
                simplify(propertyPrefactors.pre_u1 * A),
                simplify(-propertyPrefactors.pre_u2 * A),
                simplify(propertyPrefactors.pre_u3 * a0)
            )
        )
    end

    PrecomputedQuantities{Float64, Float64}(uPrefactors, ϕ, argTerms)
end

function computeUQuantitiesSym(
        filament::AFilament{0, M} where {M},
        activationsFourier::Vector{ActivationFourier{T}} where {T},
        prefactors::Vector{Prefactors}
)
    uPrefactors = Vector{Prefactors}()
    type = typeof(activationsFourier[1]).parameters[1]
    ϕ = Vector{type}()
    argTerms = Vector{type}()
    for (activationFourier, propertyPrefactors, ring) in zip(
        activationsFourier, prefactors, filament.rings)
        a0 = activationFourier.a0
        a1 = activationFourier.a1
        b1 = activationFourier.b1
        A = sqrt(a1^2 + b1^2)
        R2 = ring.geometry.R2
        α2 = ring.fiberArchitecture.α2
        push!(ϕ, simplify(-atan(b1, a1)))
        push!(argTerms, simplify(tan(α2) / R2))

        push!(
            uPrefactors,
            Prefactors(
                simplify(propertyPrefactors.pre_ζ * a0),
                simplify(propertyPrefactors.pre_u1 * A),
                simplify(-propertyPrefactors.pre_u2 * A),
                simplify(propertyPrefactors.pre_u3 * a0)
            )
        )
    end

    PrecomputedQuantities{type}(uPrefactors, ϕ, argTerms)
end

"""
    $(TYPEDSIGNATURES)

Computes the extension `ζ_hat`, and curvatures `u1_hat`,
`u2_hat`, `u3_hat` at `Z` given some `PrecomputedQuantities`.

Outputs a 4-element static `SVector` of either symbolic or numerical
(`ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`).

Used for either symbolic or numerical computation.
"""
function computeUHatSym(Z, precomputedQuantities::PrecomputedQuantities)
    ζ_hat = 1.0
    u1_hat = 0.0
    u2_hat = 0.0
    u3_hat = 0.0

    uPrefactors = precomputedQuantities.uPrefactors
    ϕ = precomputedQuantities.ϕ
    tanFactors = precomputedQuantities.argTerms
    for (uPrefactors_i, ϕ_i, tanFactor_i) in zip(uPrefactors, ϕ, tanFactors)
        ζ_hat += uPrefactors_i.pre_ζ
        u1_hat += uPrefactors_i.pre_u1 * sin(ϕ_i - Z * tanFactor_i)
        u2_hat += uPrefactors_i.pre_u2 * cos(ϕ_i - Z * tanFactor_i)
        u3_hat += uPrefactors_i.pre_u3
    end
    SVector{4}(simplify(ζ_hat), simplify(u1_hat), simplify(u2_hat), simplify(u3_hat))
end

function computeUHatSym(Z, precomputedQuantities::PrecomputedQuantitiesTapered)
    ζ_hat = 1.0
    u1_hat = 0.0
    u2_hat = 0.0
    u3_hat = 0.0

    uPrefactors = precomputedQuantities.uPrefactors
    ϕ = precomputedQuantities.ϕ
    theta2T = precomputedQuantities.argTerms
    for (uPrefactors_i, ϕ_i, theta2T_i) in zip(uPrefactors, ϕ, theta2T)
        ζ_hat += uPrefactors_i.pre_ζ(Z)
        u1_hat += uPrefactors_i.pre_u1(Z) * sin(ϕ_i - theta2T_i(Z))
        u2_hat += uPrefactors_i.pre_u2(Z) * cos(ϕ_i - theta2T_i(Z))
        u3_hat += uPrefactors_i.pre_u3(Z)
    end
    SVector{4}(ζ_hat, u1_hat, u2_hat, u3_hat)
end

"""
    $(TYPEDSIGNATURES)

Computes the extension `ζ_hat`, and curvatures `u1_hat`,
`u2_hat`, `u3_hat` at `Z` given a collection `precomp`
of precomputed quantities.

Outputs a 4-element static `SVector` of numerical
(`ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`).

Used for non-symbolic computation only.
"""
function computeUHat(Z, precomp::SMatrix{T, 6, Float64} where {T})
    ζ_hat = 1.0
    u1_hat = 0.0
    u2_hat = 0.0
    u3_hat = 0.0
    M, _ = size(precomp)

    for i in 1:M
        ζ_hat += precomp[i, 1]
        u1_hat += precomp[i, 2] * sin(precomp[i, 5] - Z * precomp[i, 6])
        u2_hat += precomp[i, 3] * cos(precomp[i, 5] - Z * precomp[i, 6])
        u3_hat += precomp[i, 4]
    end
    SVector{4}(ζ_hat, u1_hat, u2_hat, u3_hat)
end

function computeUHat(Z, precomp::SVector{M, Tuple} where {M})
    ζ_hat = 1.0
    u1_hat = 0.0
    u2_hat = 0.0
    u3_hat = 0.0
    M = length(precomp)

    for i in 1:M
        ζ_hat = ζ_hat + precomp[i][1](Z)
        u1_hat += precomp[i][2](Z) * sin(precomp[i][5] - precomp[i][6](Z))
        u2_hat += precomp[i][3](Z) * cos(precomp[i][5] - precomp[i][6](Z))
        u3_hat += precomp[i][4](Z)
    end
    SVector{4}(ζ_hat, u1_hat, u2_hat, u3_hat)
end

function computeUHat(Z, precomp::Tuple)
    ζ_hat = 1.0
    u1_hat = 0.0
    u2_hat = 0.0
    u3_hat = 0.0
    M = length(precomp)

    for i in 1:M
        ζ_hat = ζ_hat + precomp[i][1](Z)
        u1_hat += precomp[i][2](Z) * sin(precomp[i][5] - precomp[i][6](Z))
        u2_hat += precomp[i][3](Z) * cos(precomp[i][5] - precomp[i][6](Z))
        u3_hat += precomp[i][4](Z)
    end
    SVector{4}(ζ_hat, u1_hat, u2_hat, u3_hat)
end

"""
    $(TYPEDSIGNATURES)

Computes the extension `ζ_hat`, and curvatures `u1_hat`,
`u2_hat`, `u3_hat` at `Z` given a static `SMatrix` `precomp`
of precomputed quantities.

Outputs a 4-element static `SVector` of numerical
(`ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`).

Used for non-symbolic computation on the GPU.
"""
function computeUHatGPU(Z::Float32, precomp::SMatrix)
    ζ_hat = 1.0f0
    u1_hat = 0.0f0
    u2_hat = 0.0f0
    u3_hat = 0.0f0
    M, _ = size(precomp)

    for i in 1:M
        ζ_hat += precomp[i, 1]
        u1_hat += precomp[i, 2] * sin(precomp[i, 5] - Z * precomp[i, 6])
        u2_hat += precomp[i, 3] * cos(precomp[i, 5] - Z * precomp[i, 6])
        u3_hat += precomp[i, 4]
    end
    SVector{4}(ζ_hat, u1_hat, u2_hat, u3_hat)
end
#endregion ===========================
