######### Auxiliary functions for building fast functions from symbolic expressions
#region ===========================
"""
    $(TYPEDSIGNATURES)

Builds fast runtime generated functions for `ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`
based on their respective symbolic expressions.

The output also includes the extension antiderivative `ζ_hat_AD` used in computations 
with external forces.
"""
function buildUFunctions(filament::AFilament, precomputedQuantities::PrecomputedQuantities; worldage = true)
    u = computeUHat2(filament.Z, precomputedQuantities);

    if (u[1] isa PiecewiseStructure)
        uf = Vector{PiecewiseStructure}();
    else
        uf = Vector{Function}();
    end
    for ui in u
        if ui isa PiecewiseStructure
            push!(uf, eval_build_functions(ui, filament.Z, worldage = worldage));
        else
            if worldage
                push!(uf, eval(build_function(ui, filament.Z, expression = Val{false})));
            else
                push!(uf, eval(build_function(ui, filament.Z)));
            end
        end
    end

    return Tuple(uf);
end

"""
    $(TYPEDSIGNATURES)

Builds fast runtime generated functions for `ζ_hat`, `u1_hat`, `u2_hat`, `u3_hat`
based on their respective symbolic expressions.
"""
function buildUFunctionsIntrinsic(filament::AFilament, precomputedQuantities::PrecomputedQuantities; worldage = true)
    u = computeUHat2(filament.Z, precomputedQuantities);
    uf = Vector{Any}();
    for ui in u
        if ui isa PiecewiseStructure
            push!(uf, eval_build_functions(ui, filament.Z, worldage = worldage));
        else
            if worldage
                push!(uf, eval(build_function(ui, filament.Z, expression = Val{false})));
            else
                push!(uf, eval(build_function(ui, filament.Z)));
            end
        end
    end

    return Tuple(uf);
end
#endregion ===========================