######### Structures for formulaic precomputation and faster simulation
#region ===========================
"""
    $(TYPEDEF)

Numerical prefactors for the extension and curvature expressions for a given ring (for GPU computation).

$(TYPEDFIELDS)
"""
struct PrefactorsGPU
    "Extension prefactor"
    pre_ζ::Float32
    "Curvature 1 prefactor"
    pre_u1::Float32
    "Curvature 2 prefactor"
    pre_u2::Float32
    "Twist prefactor"
    pre_u3::Float32
end

"""
    $(TYPEDEF)

Symbolic or numeric prefactors for the extenion 
and curvature expressions for a given ring (symbolic or numeric computation)..

$(TYPEDFIELDS)
"""
struct Prefactors{T}
    "Extension prefactor"
    pre_ζ::T
    "Curvature 1 prefactor"
    pre_u1::T
    "Curvature 2 prefactor"
    pre_u2::T
    "Twist prefactor"
    pre_u3::T
end

"""
    $(TYPEDEF)

Precomputed quantity storage for all rings in a filament (for GPU computation).

$(TYPEDFIELDS)
"""
struct PrecomputedQuantitiesGPU
    "Vector of Prefactors for the extensions and curvatures"
    uPrefactors::Vector{Prefactors}
    "Vector of phases phi"
    ϕ::Vector{Float32}
    "Vector of tangent prefactors"
    tanFactors::Vector{Float32}
end

"""
    $(TYPEDEF)

Precomputed quantity storage for all rings in a filament (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
struct PrecomputedQuantities{U, T} <: AbstractPrecomputedQuantities
    "Vector of Prefactors for the extensions and curvatures"
    uPrefactors::Vector{Prefactors{T}}
    "Vector of phases phi"
    ϕ::Vector{U}
    "Vector of tangent prefactors"
    argTerms::Vector{T}
end

"""
    $(TYPEDEF)

Precomputed quantity storage for all rings in a filament (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
struct PrecomputedQuantitiesTapered <: AbstractPrecomputedQuantities
    "Vector of Prefactors for the extensions and curvatures"
    uPrefactors::Vector{Prefactors}
    "Vector of phases phi"
    ϕ::Vector
    "Vector of Theta2T interpolated functions"
    theta2T::Vector{Interpolations.Extrapolation}
end
#endregion ===========================