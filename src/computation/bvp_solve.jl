######### Functions for BVP solving
#region ===========================
"""
    $(TYPEDSIGNATURES)

Solves for the deformed configuration of a `filament` with 
some prescribed fibrillar activation under the influence of
gravity. Use with non-symbolic computation only.

Input parameters:
-   `m0` = initial guess for the support moments at `Z` = 0
-   `uInit` = initial guess for the `u` solution
-   `g_range` = range of gravitational acceleration values to step through

The gravitational acceleration stepping is necessary to ensure robust convergence
to a BVP solution. Typically, as little as 4 steps might be sufficient.
"""
function selfWeightSolve(filament::AFilament{T, M, A} where {T, M, A}, activation::Vector{ActivationPiecewiseGamma}, m0::Vector{Float64}, uInit::Vector{Float64}, g_range::StepRangeLen)
    activationFourier = [piecewiseGammaToFourier(activation_i) for activation_i in activation];
    selfWeightSolve(filament, activationFourier, m0, uInit, g_range)
end

function selfWeightSolve(filament::AFilament{T, M, A} where {T, M, A}, activation::Vector{ActivationFourier{T}} where T, m0::Vector{Float64}, uInit::Vector{Float64}, g_range::StepRangeLen)
    prefactors = computePropertyPrefactors(filament);
    precomputedQuantities = convertUQuantToStatic(filament, computeUQuantities(filament, activation, prefactors));

    stiffness = filament.auxiliary.stiffness;
    ρlin0Int = filament.auxiliary.ρlin0Int;
    m10::Float64, m20::Float64, m30::Float64 = m0;
    sol = 0;
    Zspan = (0.0, filament.L::Float64);
    for gi in g_range
        p = (gi, ρlin0Int, stiffness, precomputedQuantities);
        bvp = BVProblem(selfWeightDE!, selfWeightBC!, uInit, Zspan, p);

        sol = solve(bvp, Shooting(AutoVern7(Rodas4())), dt = filament.L / 100.0, abstol = 1e-12, reltol = 1e-12)
        # sol = solve(bvp, MIRK4(), dt = filament.L / 100.0, abstol = 1e-12, reltol = 1e-12);

        m10, m20, m30 = sol(0)[13:15];
        uInit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m10, m20, m30];
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

The gravitational acceleration stepping is necessary to ensure robust convergence
to a BVP solution. Typically, as little as 4 steps might be sufficient.
"""
function selfWeightSolveSym(filament::AFilament{T, M, A} where {T, M, A}, activation::Vector{ActivationFourier{T}} where T, m0::Vector{Float64}, uInit::Vector{Float64}, g_range::StepRangeLen; kwargs...)
    prefactors = computePropertyPrefactorsSym(filament);
    precomputedQuantities = computeUQuantitiesSym(filament, activation, prefactors);
    
    u_f::NTuple{4, Union{PiecewiseStructure, Function}} = buildUFunctions(filament, precomputedQuantities; kwargs...);
    m10::Float64, m20::Float64, m30::Float64 = m0;
    sol = 0;
    Zspan = (0.0, filament.L);
    for gi in g_range
        # p = (gi, filament, u_f, arcLength)
        p = (gi, filament, u_f)
        bvp = BVProblem(selfWeightDESym!, selfWeightBC!, uInit, Zspan, p)
        
        # Shooting(Tsit5()) is the default, but Vern7() has higher accuracy? Vern is definitely slightly slower though
        sol = solve(bvp, Shooting(AutoVern7(Rodas4())), dt = filament.L / 100.0, abstol = 1e-12, reltol = 1e-12) #, abstol = 1e-12, reltol = 1e-12

        m10, m20, m30 = sol(0)[13:15]
        uInit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m10, m20, m30];
    end

    sol
end

function selfWeightSolveSym(filament::AFilament{T, M, A} where {T, M, A}, activationGamma::Vector{ActivationPiecewiseGamma}, m0::Vector{Float64}, uInit::Vector{Float64}, g_range::StepRangeLen; kwargs...)
    activationFourier = Vector{ActivationFourier}([piecewiseGammaToFourier(activation) for activation in activationGamma]);
    selfWeightSolveSym(filament, activationFourier, m0, uInit, g_range; kwargs...)
end
#endregion ===========================