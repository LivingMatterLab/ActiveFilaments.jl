######### Active filament definition
#region ===========================
"""
    $(TYPEDEF)

Mechanical properties of the filament ring (for GPU computation).

$(TYPEDFIELDS)
"""
struct MechanicalPropertiesGPU <: AbstractMechanicalProperties
    "Elastic modulus"
    E::Float32
    "Poisson's ratio"
    ν::Float32
end

"""
    $(TYPEDEF)

Mechanical properties of the filament ring (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
struct MechanicalProperties <: AbstractMechanicalProperties
    "Elastic modulus"
    E::AbstractFloat
    "Poisson's ratio"
    ν::AbstractFloat
end

"""
    $(TYPEDEF)

Geometry of the filament ring (for GPU computation).

$(TYPEDFIELDS)
"""
struct GeometryGPU <: AbstractGeometry
    "Inner radius"
    R1::Float32
    "Outer radius"
    R2::Float32
end

"""
    $(TYPEDEF)

Geometry of the filament ring (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
@with_kw mutable struct Geometry{T} <: AbstractGeometry
    "Inner radius"
    R1::T
    "Outer radius"
    R2::T
    "Tapering angle at R = R2"
    phi2::Float64 = 0
end

# Convenience constructor
Geometry(R1::T, R2::T; phi2 = 0) where {T} = Geometry{T}(R1 = R1, R2 = R2, phi2 = phi2);

"""
    $(TYPEDEF)

Fiber architecture of the filament ring (for GPU computation).

$(TYPEDFIELDS)
"""
struct FiberArchitectureGPU <: AbstractFiberArchitecture
    "Helical angle"
    α2::Float32
end

"""
    $(TYPEDEF)

Fiber architecture of the filament ring (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
struct FiberArchitecture <: AbstractFiberArchitecture
    "Helical angle" # Test Union efficiency
    α2::Union{AbstractFloat, PiecewiseStructure}
end

"""
    $(TYPEDEF)

Filemant ring (for GPU computation).

$(TYPEDFIELDS)
"""
struct RingGPU <: AbstractRing
    "Ring mechanical properties"
    mechanicalProperties::MechanicalPropertiesGPU
    "Ring geometry"
    geometry::GeometryGPU
    "Fiber architecture"
    fiberArchitecture::FiberArchitectureGPU
end

"""
    $(TYPEDEF)

Filament ring (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
mutable struct Ring{T} <: AbstractRing
    "Ring mechanical properties"
    mechanicalProperties::MechanicalProperties
    "Ring geometry"
    geometry::Geometry{T}
    "Fiber architecture"
    fiberArchitecture::FiberArchitecture
end

# Convenience constructor
Ring(mechanicalProperties, geometry::Geometry{T}, fiberArchitecture) where {T} = Ring{T}(mechanicalProperties, geometry, fiberArchitecture);

"""
    $(TYPEDEF)

Inner tube of the filament (for GPU computation).

Special case of a ring with `GeometryGPU.R1 = 0` and no fibers.

$(TYPEDFIELDS)
"""
struct InnerTubeGPU <: AbstractRing
    "Inner tube mechanical properties"
    mechanicalProperties::MechanicalPropertiesGPU
    "Inner tube geometry"
    geometry::GeometryGPU
end

"""
    $(TYPEDEF)

Inner tube of the filament (symbolic or numeric computation).

Special case of a ring with `Geometry.R1 = 0` and no fibers.

$(TYPEDFIELDS)
"""
struct InnerTube{T} <: AbstractRing
    "Inner tube mechanical properties"
    mechanicalProperties::MechanicalProperties
    "Inner tube geometry"
    geometry::Geometry{T}
end

"""
    $(TYPEDEF)

Filament stiffness definition (for GPU computation).

$(TYPEDFIELDS)
"""
struct FilamentStiffnessGPU
    "Axial stiffness K0"
    K0::Float32
    "Bending stiffness K1"
    K1::Float32
    "Bending stiffness K2"
    K2::Float32
    "Torsional stiffness K3"
    K3::Float32
end

"""
    $(TYPEDSIGNATURES)

Creates a filament stiffness definition `FilamentStiffnessGPU` using an array of stiffness coefficients.
"""
FilamentStiffnessGPU(K::SVector{4, Float32}) = FilamentStiffnessGPU(K[1], K[2], K[3], K[4]);

"""
    $(TYPEDEF)

Filament stiffness definition (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
struct FilamentStiffness{T}
    "Axial stiffness K0"
    K0::T
    "Bending stiffness K1"
    K1::T
    "Bending stiffness K2"
    K2::T
    "Torsional stiffness K3"
    K3::T
end

"""
    $(TYPEDSIGNATURES)

Creates a filament stiffness definition `FilamentStiffness` using an array of stiffness coefficients.
"""
# FilamentStiffness(K::SVector{4, T}) where T = FilamentStiffness{T}(K[1], K[2], K[3], K[4]);
FilamentStiffness(K::SVector{4, Float64}) = FilamentStiffness{Float64}(K[1], K[2], K[3], K[4]);
FilamentStiffness(K::SVector{4, Vector{Float64}}) = FilamentStiffness{Vector{Float64}}(K[1], K[2], K[3], K[4]);
FilamentStiffness(K::SVector{4, T}) where T<:AbstractInterpolation = FilamentStiffness{T}(K[1], K[2], K[3], K[4]);

struct AuxiliaryProperties{A}
    ρlin0::Union{Function, Float64}
    ρlin0Int::Function
    stiffness::SVector{4, A}
end

"""
    $(TYPEDEF)

Active filament structure (for GPU computation).

Must use Float32 for all GPU computations.

$(TYPEDFIELDS)
"""
@with_kw struct AFilamentGPU
    "Length of the filament"
    L::Float32 = 1.0f0
    "Array of filament rings"
    rings::Vector{RingGPU}
    "Outer radius of the filament"
    R0::Float32 = rings[end].geometry.R2
    "Inner tube of the filament"
    innerTube::InnerTubeGPU = InnerTubeGPU(rings[1].mechanicalProperties, GeometryGPU(0, rings[1].geometry.R1))
    "Stiffness definition"
    stiffness::FilamentStiffnessGPU = FilamentStiffnessGPU(computeK(rings, innerTube))
    "Volumetric density"
    ρvol::Float32
    "Linear density"
    ρlin::Float32 = ρvol * pi * R0^2
    "Total mass"
    m::Float32 = ρvol * pi * R0^2 * L
end

"""
    $(TYPEDEF)

Active filament structure (symbolic or numeric computation).

Use Float64 if not using a GPU, since Float32 can have an impact on numerical solution accuracy.

No symbolic computation for tapered filaments is possible at this stage.

A linear tapering profile is assumed if phi2 != 0; to be generalized.

Important: Take care with very small values of phi2 (~phi2 < 0.0001), because there is a numerical
instability in Θ2T(Z, phi) for (Z, phi) -> (0, 0) since we encounter Log(1) / Sin(0) in the linear profile. 
Nonetheless, the actual configuration and the plots computed for very small phi2 will still be very close to correct 
since only the tiny portion for Z -> 0 contributes to the overall error and the applied Z-interpolations
effectively omit the instability anyway.

$(TYPEDFIELDS)
"""
@with_kw struct AFilament{V, M, A}
    "Length of the filament"
    L::Float64 = 1.0
    "Material coordinate"
    Z = LinRange(0.0, L, 32)
    "Tapering angle"
    phi2::Float64 = 0.0
    "Tapered designation"
    tapered::Bool = (phi2 != 0)
    "Array of filament rings"
    rings::Vector{Ring{T}} where T
    "Outer radius of the filament"
    R0 = rings[end].geometry.R2
    "Inner tube of the filament"
    innerTube::InnerTube = InnerTube(rings[1].mechanicalProperties, Geometry(R1 = 0.0, R2 = rings[1].geometry.R1, phi2 = phi2))
    "Stiffness definition"
    stiffness::FilamentStiffness = FilamentStiffness(computeK(rings, innerTube))
    "Volumetric density"
    ρvol::Float64
    "Total mass (if tapered, assumes linear tapering)"
    m::Float64 = tapered ? ρvol * pi / 3.0 * L * (R0[1]^2 + R0[1] * R0[end] + R0[end]^2) : ρvol * pi * R0.^2 * L
    "Auxiliary properties" # Double check if the Union here doesn't slow anything down!
    auxiliary::AuxiliaryProperties{A}
end

# Outer constructors
# function AFilament(rings::Vector{Ring{T}} where T<:AbstractFloat; L::Float64 = 1.0, ρvol::Float64 = 1000.0)
#     return AFilament{0, 1, Nothing}(L = L, rings = rings, ρvol = ρvol)
# end

# function AFilament(rings::Vector{Ring{T}} where T<:AbstractFloat; L::Float64 = 1.0, ρvol::Float64 = 1000.0, kwargs...)
#     return isempty(kwargs) ? AFilament{0, 1, Nothing}(L = L, rings = rings, ρvol = ρvol) : AFilament(rings; L = L, ρvol = ρvol, kwargs...)
# end

function AFilament(rings0::Vector{Ring{T}} where T<:AbstractFloat; 
                    N::Integer = 32, L::Float64 = 1.0, phi2::Float64 = 0.0, ρvol::Float64 = 1000.0,
                    innerTube::InnerTube = InnerTube(rings0[1].mechanicalProperties, Geometry(R1 = 0.0, R2 = rings0[1].geometry.R1, phi2 = phi2)),
                    interp::Function = cubic_spline_interpolation)
    Z = LinRange(0.0, L, N)
    @variables Zs;

    R0 = rings0[end].geometry.R2;
    if (phi2 != 0)
        f = 1 .- Z / R0 * tan(phi2);
        M = length(rings0);
        phi2_rings = Vector{Float64}(undef, M);
        rings = Vector{Ring{typeof(Z)}}(undef, M);
        for i in 1:M
            phi2_ring = atan(rings0[i].geometry.R2 / R0 * tan(phi2));
            rings[i] = Ring(rings0[i].mechanicalProperties, 
                Geometry(rings0[i].geometry.R1 * f, rings0[i].geometry.R2 * f; phi2 = phi2_ring), 
                rings0[i].fiberArchitecture);
            phi2_rings[i] = phi2_ring;
        end
        innerTube = InnerTube(innerTube.mechanicalProperties, Geometry(LinRange(0.0, 0.0, N), rings[1].geometry.R1; phi2 = atan(rings0[1].geometry.R1 / R0 * tan(phi2))))
        
        R0 = rings[end].geometry.R2;

        # expression = Val{false} prevents worldage issues when exporting with JLD2 (no "expression = Val{false}" causes the runtime generated function to not be
        # save properly in JLD2). BUT it's possible that this needs to be changed in some cases
        ρlin0 = eval(build_function(simplify(pi * ρvol * (R0[1] - Zs * tan(phi2))^2), Zs, expression = Val{false}));
        ρlin0Int = eval(build_function(simplify(pi / (3.0 * L^2) * ρvol * ((L - Zs)^3 * R0[1]^2 + (L - Zs)^2 * (L + 2 * Zs) * R0[1] * R0[end] + (L^3 - Zs^3) * R0[end]^2)), Zs, expression = Val{false}))
    else
        rings = rings0;

        ρlin0 = pi * ρvol * R0^2;

        # see other comment about worldage
        ρlin0Int = eval(build_function(ρlin0 * (L - Zs), Zs, expression = Val{false}));
    end

    stiffness = FilamentStiffness(computeK(rings, innerTube))
    stiffness_aux = phi2 != 0 ? SVector{4}([interp(Z, stiffness.K0), interp(Z, stiffness.K1), interp(Z, stiffness.K2), interp(Z, stiffness.K3)]) : 
                                SVector{4, Float64}([stiffness.K0, stiffness.K1, stiffness.K2, stiffness.K3]);
    auxProp = AuxiliaryProperties(ρlin0, ρlin0Int, stiffness_aux)
    AFilament{phi2 != 0 ? 1 : 0, length(rings), typeof(auxProp).parameters[1]}(Z = Z, L = L, rings = rings, innerTube = innerTube, ρvol = ρvol, phi2 = phi2, auxiliary = auxProp);
end

function AFilament(rings::Vector{Ring{T}} where T<:AbstractVector; 
                    N::Integer = 32, L::Float64 = 1.0, phi2 = 0.0, ρvol = 1000.0, 
                    innerTube::InnerTube = InnerTube(rings0[1].mechanicalProperties, Geometry(R1 = 0.0, R2 = rings0[1].geometry.R1, phi2 = phi2)),
                    interp::Function = cubic_spline_interpolation)
    Z = LinRange(0.0, L, N);
    @variables Zs;

    R0 = rings0[end].geometry.R2;
    if (phi2 != 0)
        # see other comment about worldage
        ρlin0 = eval(build_function(simplify(pi * ρvol * (R0[1] - Zs * tan(phi2))^2), Zs, expression = Val{false}));
        ρlin0Int = eval(build_function(simplify(pi / (3.0 * L^2) * ρvol * ((L - Zs)^3 * R0[1]^2 + (L - Zs)^2 * (L + 2 * Zs) * R0[1] * R0[end] + (L^3 - Zs^3) * R0[end]^2)), Zs, expression = Val{false}))
    else
        # see other comment about worldage
        ρlin0 = pi * ρvol * R0[1]^2;
        ρlin0Int = eval(build_function(ρlin0 * (L - Zs), Zs, expression = Val{false}));
    end

    stiffness = FilamentStiffness(computeK(rings, innerTube))
    stiffness_aux = phi2 != 0 ? SVector{4}([interp(Z, stiffness.K0), interp(Z, stiffness.K1), interp(Z, stiffness.K2), interp(Z, stiffness.K3)]) : 
                                SVector{4, Float64}([stiffness.K0, stiffness.K1, stiffness.K2, stiffness.K3]);
    auxProp = AuxiliaryProperties(ρlin0, ρlin0Int, stiffness_aux)
    AFilament{phi2 != 0 ? 1 : 0, length(rings), typeof(auxProp).parameters[1]}(Z = Z, L = L, rings = rings, innerTube = innerTube, ρvol = ρvol, phi2 = phi2, auxiliary = auxProp);
end

# function AFilament(N, L, rings0::Vector{Ring{T}} where T<:AbstractFloat, innerTube, ρvol; phi2 = 0, interp = cubic_spline_interpolation)
#     Z = LinRange(0.0, L, N)
#     @variables Zs;

#     R0 = rings0[end].geometry.R2;
#     if (phi2 != 0)
#         f = 1 .- Z / R0 * tan(phi2);
#         M = length(rings0);
#         phi2_rings = Vector{Float64}(undef, M);
#         rings = Vector{Ring{typeof(Z)}}(undef, M);
#         for i in 1:M
#             phi2_ring = atan(rings0[i].geometry.R2 / R0 * tan(phi2));
#             rings[i] = Ring(rings0[i].mechanicalProperties, 
#                 Geometry(rings0[i].geometry.R1 * f, rings0[i].geometry.R2 * f; phi2 = phi2_ring), 
#                 rings0[i].fiberArchitecture);
#             phi2_rings[i] = phi2_ring;
#         end
#         innerTube = InnerTube(innerTube.mechanicalProperties, Geometry(LinRange(0.0, 0.0, N), rings[1].geometry.R1; phi2 = atan(rings0[1].geometry.R1 / R0 * tan(phi2))))
        
#         R0 = rings[end].geometry.R2;
#         ρlin0 = eval(build_function(simplify(pi * ρvol * (R0[1] - Zs * tan(phi2))^2), Zs));
#         ρlin0Int = eval(build_function(simplify(pi / (3.0 * L^2) * ρvol * ((L - Zs)^3 * R0[1]^2 + (L - Zs)^2 * (L + 2 * Zs) * R0[1] * R0[end] + (L^3 - Zs^3) * R0[end]^2)), Zs))
#     else
#         rings = rings0;

#         ρlin0 = pi * ρvol * R0^2;
#         ρlin0Int = eval(build_function(ρlin0 * (L - Zs), Zs));
#     end

#     stiffness = FilamentStiffness(computeK(rings, innerTube))
#     # stiffness_aux = phi2 != 0 ? FilamentStiffness([interp(Z, stiffness.K0), interp(Z, stiffness.K1), interp(Z, stiffness.K2), interp(Z, stiffness.K3)]) : stiffness;
#     stiffness_aux = phi2 != 0 ? SVector{4}([interp(Z, stiffness.K0), interp(Z, stiffness.K1), interp(Z, stiffness.K2), interp(Z, stiffness.K3)]) : 
#                                 SVector{4, Float64}([stiffness.K0, stiffness.K1, stiffness.K2, stiffness.K3]);
#     auxProp = AuxiliaryProperties(ρlin0, ρlin0Int, stiffness_aux)
#     AFilament{phi2 != 0 ? 1 : 0, length(rings), typeof(auxProp).parameters[1]}(Z = Z, L = L, rings = rings, innerTube = innerTube, ρvol = ρvol, phi2 = phi2, auxiliary = auxProp);
# end

# function AFilament(N, L, rings::Vector{Ring{T}} where T<:AbstractVector, innerTube, ρvol; phi2 = 0)
#     Z = LinRange(0.0, L, N);
#     @variables Zs;

#     R0 = rings0[end].geometry.R2;
#     if (phi2 != 0)
#         ρlin0 = eval(build_function(simplify(pi * ρvol * (R0[1] - Zs * tan(phi2))^2), Zs));
#         ρlin0Int = eval(build_function(simplify(pi / (3.0 * L^2) * ρvol * ((L - Zs)^3 * R0[1]^2 + (L - Zs)^2 * (L + 2 * Zs) * R0[1] * R0[end] + (L^3 - Zs^3) * R0[end]^2)), Zs))
#     else
#         ρlin0 = pi * ρvol * R0[1]^2;
#         ρlin0Int = eval(build_function(ρlin0 * (L - Zs), Zs));
#     end

#     stiffness = FilamentStiffness(computeK(rings, innerTube))
#     # stiffness_aux = phi2 != 0 ? FilamentStiffness([interp(Z, stiffness.K0), interp(Z, stiffness.K1), interp(Z, stiffness.K2), interp(Z, stiffness.K3)]) : stiffness;
#     stiffness_aux = phi2 != 0 ? SVector{4}([interp(Z, stiffness.K0), interp(Z, stiffness.K1), interp(Z, stiffness.K2), interp(Z, stiffness.K3)]) : 
#                                 SVector{4, Float64}([stiffness.K0, stiffness.K1, stiffness.K2, stiffness.K3]);
#     auxProp = AuxiliaryProperties(ρlin0, ρlin0Int, stiffness_aux)
#     AFilament{phi2 != 0 ? 1 : 0, length(rings), typeof(auxProp).parameters[1]}(Z = Z, L = L, rings = rings, innerTube = innerTube, ρvol = ρvol, phi2 = phi2, auxiliary = auxProp);
# end
#endregion ===========================