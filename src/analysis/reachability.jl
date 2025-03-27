##########################################
### Computation of reachability clouds ###
##########################################

######### Intrinsic reachability clouds, without external loading
#region ===========================
"""
    $(TYPEDSIGNATURES)

Generates an intrinsic reachability cloud of a given `filament`
with an activation form `activationGamma`. The reachability cloud
consists of `nTrajectories` activated configurations 
computed without external loading.

The cloud data is saved to files in `path`. Setting `save_gamma_structs = false`
can be used to omit the generation of potentially large data files containing
the full gamma structures for all configurations.

"""
function generateIntrinsicReachVol(filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        nTrajectories::Int,
        path::String;
        save_gamma_structs = true)
    prefactors = computePropertyPrefactors(filament)
    M = typeof(filament).parameters[2]

    (activationsFourier, activationsGamma) = generateRandomActivations(
        activationGamma, gammaBounds, M, nTrajectories
        )
        
    println("Precomputation started...")
    @time precomputedQuantities = generatePrecomputedQuantitiesSA(
        filament,
        activationsFourier,
        prefactors,
        nTrajectories
    )
    println("Precomputation complete.")

    u0 = SVector{12, Float64}([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
    Zspan = (0.0, filament.L)
    prob = ODEProblem(
        intrinsicConfDESA,
        u0,
        Zspan,
        precomputedQuantities[1],
        save_everystep = false
    )

    prob_func = let p = precomputedQuantities
        (prob, i, repeat) -> begin
            remake(prob, p = p[i])
        end
    end

    output_func(sol, i) = (sol[end][1:3], false)

    ensemble_prob = EnsembleProblem(
        prob,
        prob_func = prob_func,
        output_func = output_func,
        safetycopy = false
    )

    println("Computing cloud...")
    @time sol = solve(
        ensemble_prob,
        Tsit5(),
        EnsembleThreads(),
        trajectories = nTrajectories,
        adaptive = false,
        dt = filament.L / 100.0
    )
    println("Cloud computed!")

    if (path != "")
        println("Saving...")
        @save string(path, ".jld2") sol

        if save_gamma_structs
            @save string(path, "_gamma.jld2") activationsGamma
        else
            n_act = [activation.N for activation in activationsGamma[1]]
            n_act_total = sum(n_act)
            running_totals = cumsum(n_act)
            idxs = [i == 1 ? (1:running_totals[1]) :
                    ((running_totals[i - 1] + 1):running_totals[i])
                    for
                    i in eachindex(running_totals)]
            activations_unrolled = MMatrix{nTrajectories, n_act_total, Float64}(undef)
            for i in 1:nTrajectories
                for j in eachindex(idxs)
                    activations_unrolled[i, idxs[j]] .= activationsGamma[i][j].γ
                end
            end
            # Stack overflow for large clouds when saving into JLD2.
            # JLD2 is not needed for raw gammas anyway.
            # @save string(path, "_gamma.jld2") activations_unrolled; 

            CSV.write(
                string(path, "_gamma.csv"),
                Tables.table(activations_unrolled),
                writeheader = false
            )
        end

        @save string(path, "_gammaBounds.jld2") gammaBounds

        x = getindex.(sol.u, 1)
        y = getindex.(sol.u, 2)
        z = getindex.(sol.u, 3)

        points = hcat(x, y, z)
        CSV.write(string(path, ".csv"), Tables.table(points), writeheader = false)
        println("Saved!")
    end

    return (sol, activationsGamma)
end

function generateIntrinsicReachVolSym(filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        nTrajectories::Int,
        path::String;
        worldage = true)
    prefactors = computePropertyPrefactors(filament)
    precomputedQuantities = generatePrecomputedQuantities(
        filament,
        activationGamma,
        gammaBounds,
        nTrajectories,
        prefactors
    )

    u_f = generateUFunctionsIntrinsic(filament, precomputedQuantities, worldage = worldage)
    println("Generated Runtime Functions")

    u0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    Zspan = (0.0, filament.L)
    prob = ODEProblem(intrinsicConfDESym!, u0, Zspan, u_f[1], save_everystep = false)

    prob_func = let p = u_f
        (prob, i, repeat) -> begin
            remake(prob, p = p[i])
        end
    end

    output_func(sol, i) = (sol[end][1:3], false)

    ensemble_prob = EnsembleProblem(
        prob,
        prob_func = prob_func,
        output_func = output_func,
        safetycopy = false
    )

    @time sol = solve(
        ensemble_prob,
        Tsit5(),
        EnsembleThreads(),
        trajectories = nTrajectories,
        adaptive = false,
        dt = filament.L / 100
    )

    @save string(path, ".jld2") sol

    @save string(path, "_gamma.jld2") activationsGammaSym

    @save string(path, "_gammaBounds.jld2") gammaBounds

    x = getindex.(sol.u, 1)
    y = getindex.(sol.u, 2)
    z = getindex.(sol.u, 3)

    points = hcat(x, y, z)
    CSV.write(string(path, ".csv"), Tables.table(points), writeheader = false)

    return sol
end
#endregion ===========================

######### Reachability clouds with external loading (self-weight)
#region ===========================
"""
    $(TYPEDSIGNATURES)

Generates a reachability cloud of a given `filament`
with an activation form `activationGamma`. The reachability cloud
consists of `nTrajectories` deformed configurations 
computed with the consideration of the filament's weight. 

Setting `save_full = true` outputs the entire BVP solution for
each configuration.

"""
function generateSelfWeightReachVol(filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        u0::Vector{Vector{Float64}},
        bcs::SVector{12, Float64},
        g_range::StepRangeLen,
        nTrajectories::Int,
        path::String;
        save_gamma_structs = true,
        save_full = false)
    prefactors = computePropertyPrefactors(filament)
    M = typeof(filament).parameters[2]

    (activationsFourier, activationsGamma) = generateRandomActivations(
        activationGamma, gammaBounds, M, nTrajectories
        )
        
    println("Precomputation started...")
    @time precomputedQuantities = generatePrecomputedQuantitiesSA(
        filament,
        activationsFourier,
        prefactors,
        nTrajectories
    )
    println("Precomputation complete.")
    stiffness = filament.auxiliary.stiffness
    ρlin0Int = filament.auxiliary.ρlin0Int

    println("Computing cloud...")
    sol = []
    for gi in g_range
        pAll = [(gi, ρlin0Int, stiffness, precomputedQuantities[i], bcs)
                for i in 1:nTrajectories]
        Zspan = (0.0, filament.L)

        bvp = TwoPointBVProblem(
                selfWeightDE!, 
                (selfWeightBCStart!, selfWeightBCEnd!), 
                u0[1], Zspan, pAll[1];
                bcresid_prototype = (zeros(12), zeros(3)))

        prob_func = let u0 = u0, p = pAll
            (prob, i, repeat) -> begin
                remake(prob, u0 = u0[i], p = p[i])
            end
        end

        if save_full
            ensemble_prob = EnsembleProblem(bvp, prob_func = prob_func, safetycopy = false)
        else
            output_func(sol, i) = ([sol[end][1:3], sol[1][13:15]], false)
            ensemble_prob = EnsembleProblem(
                bvp,
                prob_func = prob_func,
                output_func = output_func,
                safetycopy = false
            )
        end

        println(string("Started solving the cloud for g = ", gi, "..."))
        sol = solve(
            ensemble_prob,
            Shooting(AutoVern7(Rodas4())),
            EnsembleThreads(),
            trajectories = nTrajectories,
            dt = filament.L / 100.0,
            abstol = 1e-12,
            reltol = 1e-12
        )
        println("Finished solving.")

        for i in eachindex(sol)
            if save_full
                u0[i][13:15] = sol[i](0.0)[13:15]
            else
                u0[i][13:15] = [sol[i][2][1], sol[i][2][2], sol[i][2][3]]
            end
        end
        println("Done building u0")
    end
    println("Cloud computed!")

    if (path != "")
        println("Saving...")
        @save string(path, ".jld2") sol

        if save_gamma_structs
            @save string(path, "_gamma.jld2") activationsGamma
        else
            n_act = [activation.N for activation in activationsGamma[1]]
            n_act_total = sum(n_act)
            running_totals = cumsum(n_act)
            idxs = [i == 1 ? (1:running_totals[1]) :
                    ((running_totals[i - 1] + 1):running_totals[i])
                    for
                    i in eachindex(running_totals)]
            activations_unrolled = MMatrix{nTrajectories, n_act_total, Float64}(undef)
            for i in 1:nTrajectories
                for j in eachindex(idxs)
                    activations_unrolled[i, idxs[j]] .= activationsGamma[i][j].γ
                end
            end

            CSV.write(
                string(path, "_gamma.csv"),
                Tables.table(activations_unrolled),
                writeheader = false
            )
        end

        @save string(path, "_gammaBounds.jld2") gammaBounds
    end

    (sol, activationsGamma)
end

"""
    $(TYPEDSIGNATURES)

Generates a reachability cloud of a given `filament`
with an activation form `activationGamma`. The reachability cloud
consists of `nTrajectories` deformed configurations 
computed with the consideration of the filament's weight. This
function is intended for symbolic computation.

Setting `save_full = true` outputs the entire BVP solution for
each configuration.

"""
function generateSelfWeightReachVolSym(filament::AFilament,
        activationGamma::Vector{ActivationPiecewiseGamma},
        gammaBounds,
        uInit::Vector{Vector{Float64}},
        g_range::StepRangeLen,
        nTrajectories::Int;
        save_full = false)
    prefactors = computePropertyPrefactorsSym(filament)
    (precomputedQuantities, activationsGamma) = generatePrecomputedQuantitiesSym(
        filament,
        activationGamma,
        gammaBounds,
        nTrajectories,
        prefactors
    )

    u_f = generateUFunctions(filament, precomputedQuantities)
    println("Generated Runtime Functions")

    sol = []
    for gi in g_range
        pAll = [(gi, filament, u_f[i]) for i in 1:nTrajectories]
        Zspan = (0.0, filament.L)

        bvp = BVProblem(selfWeightDESym!, selfWeightBC!, uInit[1], Zspan, pAll[1])

        prob_func = let u0 = uInit, p = pAll
            (prob, i, repeat) -> begin
                remake(prob, u0 = u0[i], p = p[i])
            end
        end

        if save_full
            ensemble_prob = EnsembleProblem(bvp, prob_func = prob_func, safetycopy = false)
        else
            output_func(sol, i) = ([sol[end][1:3], sol[1][13:15]], false)
            ensemble_prob = EnsembleProblem(
                bvp,
                prob_func = prob_func,
                output_func = output_func,
                safetycopy = false
            )
        end

        println(string("Started solving for g = ", gi, "..."))
        sol = solve(
            ensemble_prob,
            Shooting(AutoVern7(Rodas4())),
            EnsembleThreads(),
            trajectories = nTrajectories,
            dt = filament.L / 20.0,
            abstol = 1e-12,
            reltol = 1e-12
        )
        println("Finished solving...")

        for i in eachindex(sol)
            if save_full
                uInit[i][13:15] = sol[i](0.0)[13:15]
            else
                uInit[i][13:15] = [sol[i][2][1], sol[i][2][2], sol[i][2][3]]
            end
        end
        println("Done building uInit")
    end

    (sol, activationsGamma)
end
#endregion ===========================