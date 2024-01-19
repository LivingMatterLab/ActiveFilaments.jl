######### Legacy types and functions
#region ===========================
struct MechanicalPropertiesSym
    E
    ν
end

struct GeometrySym
    R1
    R2
end

@with_kw struct FiberArchitectureSym
    α2
end

struct RingSym
    mechanicalProperties::MechanicalPropertiesSym
    geometry::GeometrySym
    fiberArchitecture::FiberArchitectureSym
end

struct InnerTubeSym
    mechanicalProperties::MechanicalPropertiesSym
    geometry::GeometrySym
end

struct FilamentStiffnessSym
    K0
    K1
    K2
    K3
end

FilamentStiffnessSym(K) = FilamentStiffnessSym(K[1], K[2], K[3], K[4]);

@with_kw struct AFilamentSym
    Z
    L = 1.0
    rings::Vector{RingSym}
    R0 = rings[end].geometry.R2
    innerTube::InnerTubeSym = InnerTubeSym(rings[1].mechanicalProperties, Geometry(0.0, rings[1].geometry.R1))
    stiffness::FilamentStiffnessSym = FilamentStiffnessSym(computeK(rings, innerTube))
    ρvol
    m = ρvol * pi * R0^2 * L # Change to an appropriate integral for varying R0
end

struct PrefactorsSym
    pre_ζ
    pre_u1
    pre_u2
    pre_u3
end

struct PrecomputedQuantitiesSym
    uPrefactors::Vector{PrefactorsSym}
    ϕ::Vector
    tanFactors::Vector
end

mutable struct ActivationPiecewiseGammaSym
    N::Int
    γ
    σ
    θ0
end

struct ActivationFourierSym
    a0
    a1
    b1
end

function selfWeightBCSym!(residual, u, p, t)
    residual[1] = u[1][1] - 0.0;
    residual[2] = u[1][2] - 0.0;
    residual[3] = u[1][3] - 0.0;
    residual[4] = u[1][4] - 1.0;
    residual[5] = u[1][5] - 0.0;
    residual[6] = u[1][6] - 0.0;
    residual[7] = u[1][7] - 0.0;
    residual[8] = u[1][8] - 1.0;
    residual[9] = u[1][9] - 0.0;
    residual[10] = u[1][10] - 0.0;
    residual[11] = u[1][11] - 0.0;
    residual[12] = u[1][12] - 1.0;
    residual[13] = u[end][13] - 0.0;
    residual[14] = u[end][14] - 0.0;
    residual[15] = u[end][15] - 0.0;
end
#endregion ===========================