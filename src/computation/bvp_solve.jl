#################################
### Functions for BVP solving ###
#################################
#region ===========================
"""
    $(TYPEDSIGNATURES)

Solves for the deformed configuration of a `filament` with 
some prescribed fibrillar activation under the influence of
gravity. Use with non-symbolic computation only.

Input parameters:
-   `m0` = initial guess for the support moments at `Z` = 0
-   `u0` = initial guess for the `u` solution
-   `g_range` = range of gravitational acceleration values to step through

Depending on the selected BVP solver and the applied loading, 
gravitational acceleration stepping can be necessary to ensure 
robust convergence to a BVP solution. Typically, as few as 3 steps are sufficient.
"""
function selfWeightSolve(
        filament::AFilament{T, M, A} where {T, M, A},
        activation::Vector{ActivationPiecewiseGamma},
        u0::Vector{Float64},
        bcs::SVector{12, Float64},
        g_range::StepRangeLen;
        solver::Int = 1
)
    activationFourier::Vector{ActivationFourier} = [piecewiseGammaToFourier(activation_i)
                                                    for activation_i in activation]
    selfWeightSolve(filament, activationFourier, u0, bcs, g_range; solver = solver)
end

function selfWeightSolve(
        filament::AFilament{T, M, A} where {T, M, A},
        activation::Vector{ActivationFourier},
        u0::Vector{Float64},
        bcs::SVector{12, Float64},
        g_range::StepRangeLen;
        solver::Int = 1
)
    prefactors = computePropertyPrefactors(filament)
    precomputedQuantities = convertUQuantToStatic(
        filament,
        computeUQuantities(filament, activation, prefactors)
    )

    stiffness = filament.auxiliary.stiffness
    ρlin0Int = filament.auxiliary.ρlin0Int

    sol = 0
    Zspan = (0.0, filament.L::Float64)
    for gi in g_range
        p = (gi, ρlin0Int, stiffness, precomputedQuantities, bcs)
        bvp = TwoPointBVProblem(
                    selfWeightDE!, 
                    (selfWeightBCStart!, selfWeightBCEnd!), 
                    u0, Zspan, p; 
                    bcresid_prototype = (zeros(12), zeros(3))
                )

        if solver == 1
            sol = solve(
                bvp,
                Shooting(AutoVern7(Rodas4())),
                dt = filament.L / 100.0,
                abstol = 1e-12,
                reltol = 1e-12
            )

            uInit = sol.u[1]
        elseif solver == 2
            sol = solve(
                bvp,
                MIRK4(),
                dt = filament.L / 100.0,
                abstol = 1e-12,
                reltol = 1e-12
            )

            uInit = sol.u
        end
        
        
    end

    sol
end

"""
    $(TYPEDSIGNATURES)

Solves for the deformed configuration of a `filament` with 
some prescribed fibrillar activation under the influence of
gravity. Use with either symbolic or numeric computation; 
it is most often slower than selfWeightSolve for numeric computation.

Input parameters:
-   `m0` = initial guess for the support moments at `Z` = 0
-   `uInit` = initial guess for the `u` solution
-   `g_range` = range of gravitational acceleration values to step through

Depending on the selected BVP solver, the gravitational acceleration stepping 
can be necessary to ensure robust convergence to a BVP solution. 
Typically, as little as 4 steps is sufficient.
"""
function selfWeightSolveSym(
        filament::AFilament{T, M, A} where {T, M, A},
        activation::Vector{ActivationFourier{T}} where {T},
        m0::Vector{Float64},
        uInit::Vector{Float64},
        g_range::StepRangeLen;
        kwargs...
)
    prefactors = computePropertyPrefactorsSym(filament)
    precomputedQuantities = computeUQuantitiesSym(filament, activation, prefactors)

    u_f::NTuple{4, Union{PiecewiseStructure, Function}} = buildUFunctions(
        filament, precomputedQuantities; kwargs...)
    m10::Float64, m20::Float64, m30::Float64 = m0
    sol = 0
    Zspan = (0.0, filament.L)
    for gi in g_range
        p = (gi, filament, u_f)
        bvp = BVProblem(selfWeightDESym!, selfWeightBC!, uInit, Zspan, p)

        sol = solve(
            bvp,
            Shooting(AutoVern7(Rodas4())),
            dt = filament.L / 100.0,
            abstol = 1e-12,
            reltol = 1e-12
        )

        m10, m20, m30 = sol(0)[13:15]
        uInit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m10, m20, m30]
    end

    sol
end

function selfWeightSolveSym(
        filament::AFilament{T, M, A} where {T, M, A},
        activationGamma::Vector{ActivationPiecewiseGamma},
        m0::Vector{Float64},
        uInit::Vector{Float64},
        g_range::StepRangeLen;
        kwargs...
)
    activationFourier = Vector{ActivationFourier}([piecewiseGammaToFourier(activation)
                                                   for activation in activationGamma])
    selfWeightSolveSym(filament, activationFourier, m0, uInit, g_range; kwargs...)
end
#endregion ===========================