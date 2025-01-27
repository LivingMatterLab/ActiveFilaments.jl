function convertUQuantToStatic_GPU(filament::AFilament, p::PrecomputedQuantities)
    M = length(filament.rings)
    A = MMatrix{M, 6, Float32}(zeros(M, 6))
    for j in 1:M
        A[j, :] = [
            p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
            p.uPrefactors[j].pre_u3,
            p.ϕ[j], p.tanFactors[j]]
    end
    SMatrix{M, 6}(A)
end

function convertUQuantToStatic(
        filament::AFilament{1, M} where {M},
        p::PrecomputedQuantities{Float64, Interpolations.Extrapolation}
)
    # M = length(filament.rings);
    M = typeof(filament).parameters[2]
    Tuple([(p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
               p.uPrefactors[j].pre_u3,
               p.ϕ[j], p.argTerms[j]) for j in 1:M])
end

function convertUQuantToStatic(
        filament::AFilament{0, M} where {M},
        p::PrecomputedQuantities{Float64, Float64}
)
    # M = length(filament.rings);
    M = typeof(filament).parameters[2]
    A = MMatrix{M, 6, Float64}(zeros(M, 6))
    for j in 1:M
        A[j, :] = [
            p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
            p.uPrefactors[j].pre_u3,
            p.ϕ[j], p.argTerms[j]]
    end
    SMatrix{M, 6}(A)::SMatrix{M, 6, Float64}
end