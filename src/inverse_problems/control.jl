#########################################################################
### Functions used for control calculations (inverse problem solving) ###
#########################################################################

function flattenActivation(activations::AbstractArray{<:AbstractActivation}; static = false)
    n_act = [activation.N for activation in activations]
    n_act_total = sum(n_act)
    running_totals = cumsum(n_act)
    idxs = [i == 1 ? (1:running_totals[1]) :
            ((running_totals[i - 1] + 1):running_totals[i])
            for
            i in eachindex(running_totals)]
    if static
        return activations_flattened = SVector{n_act_total, Float64}([activations[i].γ[j]
                                                                      for i in eachindex(activations)
                                                                      for
                                                                      j in eachindex(activations[i].γ)])
    else
        return activations_flattened = MVector{n_act_total, Float64}([activations[i].γ[j]
                                                                      for i in eachindex(activations)
                                                                      for
                                                                      j in eachindex(activations[i].γ)])
    end
end

# The unflattening process is incompatible with many auto-diff methods apart from FiniteDiff
function unflattenActivation!(
        activation::AbstractArray{<:AbstractActivation},
        activations_flattened::AbstractArray
)
    counter = 1
    for i in eachindex(activation)
        for j in eachindex(activation[i].γ)
            activation[i].γ[j] = activations_flattened[counter]
            counter = counter + 1
        end
    end
end

function buildDistanceFunction(
        controlObjective::ConfigurationControlObjective,
        filament::AFilament,
        activation::AbstractArray{<:AbstractActivation},
        configurationSolver::Function,
        u0 = SVector{12, Float64}([
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
        Zspan = (0.0, filament.L), args...;
        kwargs...)
    if (controlObjective.propertyType === r)
        return (function f(x, p)
            unflattenActivation!(activation, x)
            sol = configurationSolver(filament, activation, u0, Zspan, args...)

            # Adapt for vectors
            sol_r = sol(controlObjective.args[1])[1:3]
            return sqeuclidean(sol_r, controlObjective.properties[1, :])
        end)
    else
        # Other property types are not used in current code. Implement as needed. 
    end
end

"""
    $(TYPEDSIGNATURES)

Optimizes the fibrillar activation in a `filament` given a `controlObjective`
and an `activation_structure`. The initial guess for the optimization is
passed as input in `activation0`.

"""
function optimizeActivation(
        controlObjective::ConfigurationControlObjective,
        filament::AFilament{V, M, A} where {V, M, A},
        activation_structure::AbstractArray{<:AbstractActivation},
        activation0::AbstractArray{<:AbstractActivation},
        args...;
        kwargs...
)
    activation = deepcopy(activation_structure)
    f = buildDistanceFunction(
        controlObjective,
        filament,
        activation,
        selfWeightSolve,
        args...
    )
    x0 = flattenActivation(activation0)
    prob = OptimizationProblem(f, x0, [1.0, 2.0]; kwargs...)
    sol = solve(prob, NLopt.G_MLSL(), local_method = NLopt.LN_NELDERMEAD(), maxtime = 10.0)

    sol
end

function optimizeActivation(
        controlObjective::ConfigurationControlObjective,
        filament::AFilament,
        activation_structure::AbstractArray{<:AbstractActivation},
        activations0::AbstractArray{<:AbstractArray{<:AbstractActivation}},
        args...;
        kwargs...
)
    f = buildDistanceFunction(
        controlObjective,
        filament,
        activation_structure,
        solveIntrinsic,
        args...
    )
    f = OptimizationFunction(f, Optimization.AutoFiniteDiff())

    sols = Vector(undef, length(activations0))
    Threads.@threads for i in eachindex(activations0)
        x0 = flattenActivation(activations0[i])
        prob = OptimizationProblem(f, x0, [1.0, 2.0]; kwargs...)
        sols[i] = solve(prob, Optim.NelderMead())
    end

    sols
end

function distance_function(
        x::Vector{Float64},
        control_objective::ConfigurationControlObjective,
        trunk::TrunkFast{T, N},
        bvp::BVProblem;
        m0::Vector{Float64} = [0.0, 0.0, 0.0],
        u0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.0, 1.0, m0[1], m0[2], m0[3]],
        Zspan::Tuple{Float64, Float64} = (0.0, trunk.trunk.L),
        abstol::Float64 = 1e-8,
        reltol::Float64 = 1e-6,
        ζ_hat_low::Float64 = 0.67,
        ζ_hat_high::Float64 = 1.35,
        solver = 1,
        kwargs...
) where {T, N}
    arg = control_objective.args
    prop = control_objective.properties
    types = control_objective.propertyTypes
    weights = control_objective.weights

    X = [x[1:12]..., 0.0, x[13:26]..., 0.0, x[27], x[28]]
    γ = (
        SMatrix{3, 5, Float64}(transpose(reshape(X[1:15], (5, 3)))),
        SMatrix{3, 5, Float64}(transpose(reshape(X[16:30], (5, 3))))
    )

    # Include rotation as DOFs                  
    bc_v = rotate_bc(trunk, (x[29], x[30]))

    new_bcs = (x[29] == bvp.p[6][1] && x[30] == bvp.p[6][2]) ? nothing :
              SVector{12, Float64}(bc_v)

    sol, a = self_weight_solve_single(bvp, trunk, γ;
        uInit = u0,
        new_bcs = new_bcs,
        abstol = abstol, reltol = reltol,
        solver = solver,
        kwargs...)

    out = 0.0
    for i in eachindex(prop)
        prop_i = prop[i]
        arg_i = arg[i]
        t = types[i]
        weights_i = weights[i]
        for j in eachindex(arg_i)
            sol_arg = sol(arg_i[j])
            if t === r
                sol_prop = sol_arg[1:3]
                dist_r = euclidean(sol_prop, prop_i[j, :])
                out = out + weights_i[j] * dist_r
            elseif t === d3
                sol_prop = sol_arg[10:12]
                dist_d3 = euclidean(sol_prop, prop_i[j, :])
                out = out + weights_i[j] * dist_d3
            end
        end
    end
    ζ_hat_penalty = 1e3 * (
        any(a.u_hat_array[4][:, 1] .< ζ_hat_low) ||
        any(a.u_hat_array[4][:, 1] .> ζ_hat_high)
    )
    out = out + ζ_hat_penalty
    @info "Current objective: $out"
    return out
end

function build_distance_function(
        control_objective::ConfigurationControlObjective,
        trunk::TrunkFast{T, N},
        configurationSolver::Function,
        bvp::BVProblem,
        x_previous::Vector{Float64},
        m0::Vector{Float64} = [0.0, 0.0, 0.0],
        uInit::Vector{Float64} = [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            m0[1],
            m0[2],
            m0[3]
        ],
        Zspan = (0.0, trunk.trunk.L), args...; optimize_bc = false,
        abstol = 1e-8,
        reltol = 1e-6,
        kwargs...) where {T, N}
    arg = control_objective.args
    prop = control_objective.properties
    types = control_objective.propertyTypes
    weights = control_objective.weights
    if optimize_bc
        return (
            (x, p) -> begin
            X = [x[1:12]..., 0.0, x[13:26]..., 0.0, x[27], x[28]]
            γ = (
                SMatrix{3, 5, Float64}(transpose(reshape(X[1:15], (5, 3)))),
                SMatrix{3, 5, Float64}(transpose(reshape(X[16:30], (5, 3))))
            )

            # Include rotation as DOFs                  
            bc_v = rotate_bc(trunk, (x[29], x[30]))

            new_bcs = (x[29] == bvp.p[6][1] && x[30] == bvp.p[6][2]) ? nothing :
                      SVector{12, Float64}(bc_v)

            sol, _ = configurationSolver(
                bvp,
                trunk,
                γ,
                args...;
                uInit = [bc_v..., 0.0, 0.0, 0.0],
                new_bcs = new_bcs,
                abstol = abstol,
                reltol = reltol,
                kwargs...
            )

            out = 0.0
            for i in eachindex(prop)
                prop_i = prop[i]
                arg_i = arg[i]
                t = types[i]
                weights_i = weights[i]
                for j in eachindex(arg_i)
                    sol_arg = sol(arg_i[j])
                    if t === r
                        sol_prop = sol_arg[1:3]
                        dist_r = euclidean(sol_prop, prop_i[j, :])
                        out = out + weights_i[j] * dist_r
                    elseif t === d3
                        sol_prop = sol_arg[10:12]
                        dist_d3 = euclidean(sol_prop, prop_i[j, :])
                        out = out + weights_i[j] * dist_d3
                    end
                end
            end
            @info "Current objective: $out"
            return out
        end
        )
    else
        return (x, p) -> begin
            X = [x[1:12]..., 0.0, x[13:26]..., 0.0, x[27], x[28]]
            γ = (
                SMatrix{3, 5, Float64}(transpose(reshape(X[1:15], (5, 3)))),
                SMatrix{3, 5, Float64}(transpose(reshape(X[16:30], (5, 3))))
            )

            sol, _ = configurationSolver(bvp, trunk, γ, args...; uInit = bvp.u0)

            out = 0.0
            for i in eachindex(prop)
                prop_i = prop[i]
                arg_i = arg[i]
                t = types[i]
                weights_i = weights[i]
                for j in eachindex(arg_i)
                    sol_arg = sol(arg_i[j])
                    if t === r
                        sol_prop = sol_arg[1:3]
                    elseif t === d3
                        sol_prop = sol_arg[10:12]
                    end

                    out = out + weights_i[j] * euclidean(sol_prop, prop_i[j, :])
                end
            end
            @info "Current objective: $out"
            return out
        end
    end
end

function build_distance_function(
        control_objective::ConfigurationControlObjective,
        trunk::TrunkFast{T, N},
        ivp::ODEProblem,
        u0::SVector{12, Float64} = (@SVector [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0
        ]),
        Zspan = (0.0, trunk.trunk.L), args...;
        ζ_hat_low = 0.67,
        ζ_hat_high = 1.35,
        kwargs...) where {T, N}
    arg = control_objective.args
    prop = control_objective.properties
    types = control_objective.propertyTypes
    weights = control_objective.weights

    return (
        (x, p) -> begin
        X = [x[1:12]..., 0.0, x[13:26]..., 0.0, x[27], x[28]]
        γ = (
            SMatrix{3, 5, Float64}(transpose(reshape(X[1:15], (5, 3)))),
            SMatrix{3, 5, Float64}(transpose(reshape(X[16:30], (5, 3))))
        )

        # Include rotation as DOFs                  
        bc_v = rotate_bc(trunk, (x[29], x[30]))

        sol, a = ivp_solve_single(ivp, trunk, γ, args...; new_u0 = bc_v, kwargs...)

        out = 0.0
        for i in eachindex(prop)
            prop_i = prop[i]
            arg_i = arg[i]
            t = types[i]
            weights_i = weights[i]
            for j in eachindex(arg_i)
                sol_arg = sol(arg_i[j])
                if t === r
                    sol_prop = sol_arg[1:3]
                    dist_r = euclidean(sol_prop, prop_i[j, :])
                    out = out + weights_i[j] * dist_r
                elseif t === d3
                    sol_prop = sol_arg[10:12]
                    dist_d3 = euclidean(sol_prop, prop_i[j, :])
                    out = out + weights_i[j] * dist_d3
                end
            end
        end
        ζ_hat_penalty = 1e3 * (
            any(a.u_hat_array[4][:, 1] .< ζ_hat_low) ||
            any(a.u_hat_array[4][:, 1] .> ζ_hat_high)
        )
        out = out + ζ_hat_penalty
        # @info "Current objective: $out"
        return out
    end
    )
end

function rotate_bc(trunk::TrunkFast{T, N}, θ::Tuple{Float64, Float64}) where {T, N}
    rot1 = AngleAxis(θ[1], 0.0, 0.0, 1.0) # Rotation around Z is most intuitive in this case
    rot2 = AngleAxis(θ[2], 0.0, 1.0, 0.0)
    rot = rot1 * rot2 # Rotation is intrinsic and not extrinsic (1 * 2 instead of 2 * 1)

    cond = trunk.trunk.clamping_condition
    sphere_joint = trunk.trunk.sphere_joint
    sphere = sphere_joint.sphere

    bc = hcat(cond.d10, cond.d20, cond.d30)

    bc = rot * bc

    r0 = sphere.c + rot * sphere_joint.sD_hat * sphere.r

    bc_v = SVector{12, Float64}(vcat(r0, reshape(bc, 9)))

    bc_v
end

function optimize_activation(
        control_objective::ConfigurationControlObjective,
        trunk::TrunkFast{T, N},
        bvp::BVProblem,
        args...;
        optimize_bc = false,
        x0::Vector{Float64} = (optimize_bc ? zeros(30) : zeros(28)),
        x_previous::Vector{Float64} = x0,
        u0 = nothing,
        maxtime = 60.0,
        abstol = 1e-8,
        reltol = 1e-6,
        local_min = false,
        solver = 1,
        kwargs...) where {T, N}
    if !isnothing(u0)
        bvp = remake(bvp; u0 = u0)
    end

    f = (x, p) -> distance_function(x, control_objective, trunk, bvp; u0 = u0,
        abstol = abstol, reltol = reltol, solver = solver, kwargs...)

    prob = OptimizationProblem(f, x0, nothing; kwargs...)

    println("Optimization start")

    if local_min
        sol = solve(prob, NLopt.LN_SBPLX(), maxtime = maxtime)
    else
        sol = solve(prob, OptimizationNOMAD.NOMADOpt(), maxtime = maxtime)
    end

    @info ":::::::::::::: FINAL OBJECTIVE = $(sol.objective) ::::::::::::::"

    X = [sol[1:12]..., 0.0, sol[13:26]..., 0.0, sol[27], sol[28]]
    γ = (
        SMatrix{3, 5, Float64}(transpose(reshape(X[1:15], (5, 3)))),
        SMatrix{3, 5, Float64}(transpose(reshape(X[16:30], (5, 3))))
    )

    if optimize_bc
        θ = (sol[29], sol[30])
        new_bcs = SVector{12, Float64}(rotate_bc(trunk, θ))
        return sol[1:28], γ, θ, new_bcs, sol.objective
    else
        return sol[1:28], γ
    end
end

function optimize_activation(
        control_objective::ConfigurationControlObjective,
        trunk::TrunkFast{T, N},
        ivp::ODEProblem, args...;
        x0::Union{SVector{30, Float64}, Vector{Float64}} = (@SVector zeros(30)),
        u0 = nothing,
        maxtime = 60.0,
        local_min = false,
        alt_method = false,
        ζ_hat_high = 1.35,
        kwargs...) where {T, N}
    if !isnothing(u0)
        ivp = remake(ivp; u0 = u0)
    end

    f = build_distance_function(
        control_objective,
        trunk,
        ivp,
        args...;
        ζ_hat_high = ζ_hat_high,
        kwargs...
    )

    println(f(x0, nothing))

    prob = OptimizationProblem(f, x0, nothing; kwargs...)

    println("Optimization start")

    if alt_method
        sol = solve(
            prob,
            OptimizationBBO.BBO_adaptive_de_rand_1_bin_radiuslimited(),
            maxtime = maxtime;
            NThreads = Threads.nthreads() - 1
        )
    else
        if (local_min)
            sol = solve(prob, NLopt.LN_SBPLX(), maxtime = maxtime)
        else
            sol = solve(
                prob,
                NLopt.G_MLSL_LDS(),
                local_method = NLopt.LN_NELDERMEAD(),
                maxtime = maxtime
            )
        end
    end

    @info "::::::::::::FINAL OBJECTIVE = $(sol.objective) ::::::::::::::"

    X = [sol[1:12]..., 0.0, sol[13:26]..., 0.0, sol[27], sol[28]]
    γ = (
        SMatrix{3, 5, Float64}(transpose(reshape(X[1:15], (5, 3)))),
        SMatrix{3, 5, Float64}(transpose(reshape(X[16:30], (5, 3))))
    )

    θ = (sol[29], sol[30])
    new_bcs = SVector{12, Float64}(rotate_bc(trunk, θ))
    return sol[1:28], γ, θ, new_bcs, sol.objective
end