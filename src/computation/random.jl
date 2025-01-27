####################################
### Random activation generators ###
####################################
"""
    $(TYPEDSIGNATURES)

Generates random activations from a uniform distribution U(`a`, `b`)^3
of the Fourier coefficients `ActivationFourier.a0`, `ActivationFourier.a1`,
`ActivationFourier.b1`.
"""
function generateRandomActivations(a::Float32, b::Float32, M::Int, n::Int)
    activationsFourier = Vector{Vector{ActivationFourier}}()
    for i in 1:n
        activationFourier = Vector{ActivationFourier}()
        for j in 1:M
            fourierCoeff = rand(Uniform(a, b), 3)
            push!(
                activationFourier,
                ActivationFourier(fourierCoeff[1], fourierCoeff[2], fourierCoeff[3])
            )
        end
        push!(activationsFourier, activationFourier)
    end
    return activationsFourier
end

"""
    $(TYPEDSIGNATURES)

Copies the `ActivationPiecewiseGamma` `struct`.
"""
Base.copy(a::ActivationPiecewiseGamma) = ActivationPiecewiseGamma(a.N, a.γ, a.σ, a.θ0);

function generateRandomActivations(
        activationGamma::ActivationPiecewiseGamma,
        a::Float32,
        b::Float32,
        M::Int,
        n::Int
)
    activationsFourier = Vector{Vector{ActivationFourier}}()
    empty!(activationsGamma)
    for i in 1:n
        activationFourier = Vector{ActivationFourier}()
        for j in 1:M
            γ = rand(Uniform(a, b), activationGamma.N)
            activationGamma.γ = γ
            push!(activationsGamma, copy(activationGamma))
            activationFourierRing = piecewiseGammaToFourier(activationGamma)
            push!(activationFourier, activationFourierRing)
        end
        push!(activationsFourier, activationFourier)
    end
    return activationsFourier
end

function generateRandomActivations(
        activationGammaParent::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        M::Int,
        n::Int;
)
    activationsFourier = Vector{Vector{ActivationFourier}}()
    activationsGamma = Vector{Vector{ActivationPiecewiseGamma}}()

    for i in 1:n
        activationFourier = Vector{ActivationFourier}()
        activationGamma = Vector{ActivationPiecewiseGamma}()
        for j in 1:M
            N = activationGammaParent[j].N
            γ = Vector{Float64}(undef, N)
            for k in 1:N
                if gammaBounds[j][1][k] != gammaBounds[j][2][k]
                    γ[k] = rand(Uniform(gammaBounds[j][1][k], gammaBounds[j][2][k]))
                else
                    γ[k] = gammaBounds[j][1][k]
                end
            end
            activationGammaRing = copy(activationGammaParent[j])
            activationGammaRing.γ = γ
            push!(activationGamma, activationGammaRing)
            activationFourierRing = piecewiseGammaToFourier(activationGammaRing)
            push!(activationFourier, activationFourierRing)
        end
        push!(activationsGamma, activationGamma)
        push!(activationsFourier, activationFourier)
    end
    return (activationsFourier, activationsGamma)
end

function generatePrecomputedQuantities(
        filament::AFilament,
        a::Float32,
        b::Float32,
        n::Int,
        prefactors::Vector{Prefactors}
)
    activationsFourier = generateRandomActivations(a, b, length(filament.rings), n)
    precomputedQuantities = Vector{PrecomputedQuantities}()
    for i in 1:n
        push!(
            precomputedQuantities,
            computeUQuantities2(filament, activationsFourier[i], prefactors)
        )
    end
    precomputedQuantities
end

function generatePrecomputedQuantities(
        filament::AFilament,
        activationGamma::ActivationPiecewiseGamma,
        a::Float32,
        b::Float32,
        n::Int,
        prefactors::Vector{Prefactors}
)
    activationsFourier = generateRandomActivations(
        activationGamma, a, b, length(filament.rings), n)
    precomputedQuantities = Vector{PrecomputedQuantities}()
    for i in 1:n
        push!(
            precomputedQuantities,
            computeUQuantities2(filament, activationsFourier[i], prefactors)
        )
    end
    precomputedQuantities
end

function generatePrecomputedQuantities(
        filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        n::Int,
        prefactors::Vector{Prefactors}
)
    M = length(filament.rings)
    (activationsFourier, activationsGamma) = generateRandomActivations(
        activationGamma, gammaBounds, M, n)
    precomputedQuantities = Vector{PrecomputedQuantities}()
    for i in 1:n
        push!(
            precomputedQuantities,
            computeUQuantities2(filament, activationsFourier[i], prefactors)
        )
    end
    return (precomputedQuantities, activationsGamma)
end

function generatePrecomputedQuantitiesSym(
        filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        n::Int,
        prefactors::Vector{Prefactors}
)
    M = length(filament.rings)
    (activationsFourier, activationsGamma) = generateRandomActivations(
        activationGamma, gammaBounds, M, n)
    precomputedQuantities = Vector{PrecomputedQuantities}()
    for i in 1:n
        push!(
            precomputedQuantities,
            computeUQuantitiesSym(filament, activationsFourier[i], prefactors)
        )
    end
    return (precomputedQuantities, activationsGamma)
end

function generatePrecomputedQuantitiesGPU(
        filament::AFilament,
        a::Float32,
        b::Float32,
        n::Int,
        prefactors::Vector{Prefactors}
)
    M = length(filament.rings)
    activationsFourier = generateRandomActivations(a, b, M, n)
    precomputedQuantities = Vector{SMatrix{M, 6, Float32}}()
    for i in 1:n
        p = computeUQuantities2(filament, activationsFourier[i], prefactors)
        A = MMatrix{M, 6, Float32}(zeros(M, 6))
        for j in 1:M
            A[j, :] = [
                p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
                p.uPrefactors[j].pre_u3,
                p.ϕ[j], p.tanFactors[j]]
        end
        push!(precomputedQuantities, SMatrix{M, 6}(A))
    end
    precomputedQuantities
end

function generatePrecomputedQuantitiesGPU(
        filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        n::Int,
        prefactors::Vector{Prefactors}
)
    M = length(filament.rings)
    activationsFourier = generateRandomActivations(activationGamma, gammaBounds, M, n)
    precomputedQuantities = Vector{SMatrix{M, 6, Float32}}()
    for i in 1:n
        p = computeUQuantities2(filament, activationsFourier[i], prefactors)
        A = MMatrix{M, 6, Float32}(zeros(M, 6))
        for j in 1:M
            A[j, :] = [
                p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
                p.uPrefactors[j].pre_u3,
                p.ϕ[j], p.tanFactors[j]]
        end
        push!(precomputedQuantities, SMatrix{M, 6}(A))
    end
    precomputedQuantities
end

"""
    $(TYPEDSIGNATURES)

Precomputes and returns all quantities necessary 
to maximize the speedup of the reachability cloud computation.

The function can take as input either tapered or non-tapered
structures and either gamma or Fourier activation forms.

If Vector{ActivationPiecewiseGamma} is passed as input, then a gammaBounds
array has to be provided as well. The function will then sample random activations
based on the form of activationGamma.

If the Fourier activation format is needed instead, the function takes as input 
Vector{Vector{ActivationFourier}} of pre-generated activations and uses
those to precompute the necessary quantities.

"""
function generatePrecomputedQuantitiesSA(
        filament::AFilament{1, M, A} where {M, A},
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        n::Int,
        prefactors::Vector{Prefactors}
)
    M = typeof(filament).parameters[2]
    (activationsFourier, activationsGamma) = generateRandomActivations(
        activationGamma, gammaBounds, M, n)

    p_all = [computeUQuantities(filament, activationsFourier[i], prefactors) for i in 1:n]

    return Tuple([Tuple([(p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1,
                             p.uPrefactors[j].pre_u2,
                             p.uPrefactors[j].pre_u3,
                             p.ϕ[j], p.argTerms[j]) for j in 1:M]) for p in p_all]),
            activationsGamma
end

function generatePrecomputedQuantitiesSA(
        filament::AFilament{0, M, A} where {M, A},
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        n::Int,
        prefactors::Vector{Prefactors}
)
    M = typeof(filament).parameters[2]
    (activationsFourier, activationsGamma) = generateRandomActivations(
        activationGamma, gammaBounds, M, n)

    precomputedQuantities = Vector{SMatrix{M, 6, Float64}}()
    for i in 1:n
        p = computeUQuantities(filament, activationsFourier[i], prefactors)
        A = MMatrix{M, 6, Float64}(zeros(M, 6))
        for j in 1:M
            A[j, :] = [
                p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
                p.uPrefactors[j].pre_u3,
                p.ϕ[j], p.argTerms[j]]
        end
        push!(precomputedQuantities, SMatrix{M, 6}(A))
    end
    (precomputedQuantities, activationsGamma)
end

function generatePrecomputedQuantitiesSA(
        filament::AFilament{1, M, A} where {M, A},
        activationsFourier::Vector{Vector{ActivationFourier}},
        prefactors::Vector{Prefactors},
        n::Int
)
    M = typeof(filament).parameters[2]
    p_all = Vector{PrecomputedQuantities{Float64, Interpolations.Extrapolation}}(undef, n)
    Threads.@threads for i in eachindex(activationsFourier)
        p_all[i] = computeUQuantities(filament, activationsFourier[i], prefactors)
    end

    return Tuple([Tuple([(p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1,
                             p.uPrefactors[j].pre_u2,
                             p.uPrefactors[j].pre_u3,
                             p.ϕ[j], p.argTerms[j]) for j in 1:M]) for p in p_all])
end

function generatePrecomputedQuantitiesSA(
        filament::AFilament{0, M, A} where {M, A},
        activationsFourier::Vector{Vector{ActivationFourier}},
        prefactors::Vector{Prefactors},
        n::Int
)
    M = typeof(filament).parameters[2]
    p_all = Vector{SMatrix{M, 6, Float64}}(undef, n)
    Threads.@threads for i in eachindex(activationsFourier)
        p = computeUQuantities(filament, activationsFourier[i], prefactors)
        A = MMatrix{M, 6, Float64}(zeros(M, 6))
        for j in 1:M
            A[j, :] = [
                p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
                p.uPrefactors[j].pre_u3,
                p.ϕ[j], p.argTerms[j]]
        end
        p_all[i] = SMatrix{M, 6}(A)
    end

    return p_all
end

function generateActivations(activationsGamma::Vector{Vector{ActivationPiecewiseGamma}}, M)
    activationsFourier = Vector{Vector{ActivationFourier}}()
    for i in 1:length(activationsGamma)
        activationFourier = Vector{ActivationFourier}()
        for j in 1:M
            activationGammaRing = activationsGamma[i][j]
            activationFourierRing = piecewiseGammaToFourier(activationGammaRing)
            push!(activationFourier, activationFourierRing)
        end
        push!(activationsFourier, activationFourier)
    end
    return activationsFourier
end

function generatePrecomputedQuantitiesAct(
        filament::AFilament,
        activationsGamma::Vector{Vector{ActivationPiecewiseGamma}},
        prefactors::Vector{Prefactors}
)
    M = length(filament.rings)
    activationsFourier = generateActivations(activationsGamma, M)
    precomputedQuantities = Vector{SMatrix{M, 6, Float32}}()
    for i in eachindex(activationsGamma)
        p = computeUQuantities2(filament, activationsFourier[i], prefactors)
        A = MMatrix{M, 6, Float32}(zeros(M, 6))
        for j in 1:M
            A[j, :] = [
                p.uPrefactors[j].pre_ζ, p.uPrefactors[j].pre_u1, p.uPrefactors[j].pre_u2,
                p.uPrefactors[j].pre_u3,
                p.ϕ[j], p.tanFactors[j]]
        end
        push!(precomputedQuantities, SMatrix{M, 6}(A))
    end
    precomputedQuantities
end

function generateUFunctions(
        filament::AFilament,
        precomputedQuantities::Vector{PrecomputedQuantities};
        worldage = true
)
    u_f = []
    for quantities in precomputedQuantities
        push!(u_f, buildUFunctions(filament, quantities, worldage = worldage))
    end
    return Tuple(u_f)
end

function generateUFunctionsIntrinsic(
        filament::AFilament,
        precomputedQuantities::Vector{PrecomputedQuantities};
        worldage = true
)
    u_f = []
    for quantities in precomputedQuantities
        push!(u_f, buildUFunctionsIntrinsic(filament, quantities, worldage = worldage))
    end
    return Tuple(u_f)
end