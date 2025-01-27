###############################################
### Functions for IVP configuration solving ###
###############################################
#region ===========================
"""
    $(TYPEDSIGNATURES)

Solves for the activated configuration of a `filament` with 
some prescribed fibrillar activation.

Input parameters:
-   `u0` = initial condition for the `u` solution
-   `Zspan` = `Tuple` indicating the first and last `Z` value

"""
function solveIntrinsic(
        filament::AFilament,
        activation::Vector{ActivationFourier},
        u0 = SVector{12, Float64}([
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
        Zspan = (0.0, filament.L);
        kwargs...
)
    prefactors = computePropertyPrefactors(filament)
    precomputedQuantities = convertUQuantToStatic(
        filament,
        computeUQuantities(filament, activation, prefactors; kwargs...)
    )
    prob = ODEProblem(intrinsicConfDESA, u0, Zspan, precomputedQuantities)

    sol = solve(
        prob,
        AutoVern7(Rodas4()),
        dt = filament.L / 100.0,
        abstol = 1e-12,
        reltol = 1e-12
    )
    sol
end

function solveIntrinsic(
        filament::AFilament,
        activation::Vector{ActivationPiecewiseGamma},
        u0 = SVector{12, Float64}([
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
        Zspan = (0.0, filament.L);
        kwargs...
)
    activationFourier = Vector{ActivationFourier}([piecewiseGammaToFourier(activation_i)
                                                   for activation_i in activation])
    solveIntrinsic(filament, activationFourier, u0, Zspan; kwargs...)
end
#endregion ===========================