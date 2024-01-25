module PlottingExt
using ActiveFilaments
using ColorSchemes
using CairoMakie
using GLMakie
using Plots
using Interpolations
using Rotations
using StaticArrays

# import ActiveFilaments: plotReachabilityCloudRGB, plotReachabilityCloud

function ActiveFilaments.plotReachabilityCloudRGB(sols, activationsGamma::Matrix{Float64}, gammaBounds, axesLimits; 
    gravity = false, full_solution = false, flipped = false, showBox = true, hideAll = false, 
    azimuth = 1.275 * pi, elevation = pi / 8, 
    resolution = (3840, 2160), 
    markersize = 2,
    perspectiveness = 0.0,
    transparency = false,
    opacity = 1.0,
    kwargs...)

    if full_solution
        if gravity
            points = [sol_i[end][1:3] for sol_i in sols]; # Assuming that second element of each sol is the support moment vector
        else
            points = sols.u; # ????
        end
    else
        if gravity
            points = [sol_i[1] for sol_i in sols]; # Assuming that second element of each sol is the support moment vector
        else
            points = sols.u;
        end
    end

    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    activations = [abs.(activationsGamma[i, :]) for i in axes(activationsGamma, 1)];
    
    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}();
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([gammaBounds[j][1][k], gammaBounds[j][2][k]]))));
            push!(gammaBoundsPairs, bounds);
        end
    end

    GLMakie.activate!()
    fig = Figure(resolution = resolution)
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = azimuth, elevation = elevation, perspectiveness = perspectiveness; kwargs...)
    GLMakie.scatter!(x, y, z, markersize = markersize, 
                    color = [
                        RGBA(activations[i][1] / gammaBoundsPairs[1][2], 
                        activations[i][2] / gammaBoundsPairs[2][2], 
                        activations[i][3] / gammaBoundsPairs[3][2], 
                        opacity) 
                        for i in eachindex(activations)
                            ], 
                        transparency = transparency)

    if (flipped)
        limits!(ax, -axesLimits[2], -axesLimits[1], axesLimits[3], axesLimits[4], -axesLimits[6], -axesLimits[5])
    else
        limits!(ax, axesLimits[1], axesLimits[2], axesLimits[3], axesLimits[4], axesLimits[5], axesLimits[6])
    end

    # if (hideAll)
    #     hidedecorations!(ax)
    #     hidespines!(ax)
    # end

    if (showBox)
        hidedecorations!(ax, grid = false);
    else
        hidedecorations!(ax);
        hidespines!(ax);
    end

    fig
end

function ActiveFilaments.plotReachabilityCloud(sols, activationsGamma::Matrix{Float64}, activationIndex::Integer, gammaBounds, axesLimits;
                                gravity = false, full_solution = false, flipped = false, showBox = true, hideAll = false,
                                azimuth = 1.275 * pi, elevation = pi / 8, 
                                resolution = (3840, 2160), 
                                markersize = 2,
                                perspectiveness = 0.0,
                                transparency = false,
                                opacity = 1.0,
                                colormap = :thermal,
                                maxAbsGamma = nothing,
                                maxAbsGammaIndexed = false,
                                kwargs...)
    if full_solution
        if gravity
            points = [sol_i[end][1:3] for sol_i in sols]; # Assuming that second element of each sol is the support moment vector
        else
            points = sols.u; # ????
        end
    else
        if gravity
            points = [sol_i[1] for sol_i in sols]; # Assuming that second element of each sol is the support moment vector
        else
            points = sols.u;
        end
    end

    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    activations = [abs.(activationsGamma[:, i]) for i in axes(activationsGamma, 2)];

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}();
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([gammaBounds[j][1][k], gammaBounds[j][2][k]]))));
            push!(gammaBoundsPairs, bounds);
        end
    end

    if !isnothing(maxAbsGamma)
        if (maxAbsGammaIndexed)
            activationMap = activations[activationIndex] .<= maxAbsGamma;
        else
            activationMap = (activations[1] .<= maxAbsGamma)
            for i in 2:size(activationsGamma, 2)
                activationMap = activationMap .& (activations[i] .<= maxAbsGamma);
            end
        end

        x = x[activationMap];
        y = y[activationMap];
        z = z[activationMap];

        for i in axes(activationsGamma, 2)
            activations[i] = activations[i][activationMap];
        end
    end

    GLMakie.activate!()
    fig = Figure(resolution = resolution)
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = azimuth, elevation = elevation, perspectiveness = perspectiveness; kwargs...)
    GLMakie.scatter!(x, y, z, markersize = markersize, color = activations[activationIndex], colormap = (colormap, opacity), colorrange = gammaBoundsPairs[activationIndex], transparency = transparency)
    
    if (flipped)
        limits!(ax, -axesLimits[2], -axesLimits[1], axesLimits[3], axesLimits[4], -axesLimits[6], -axesLimits[5])
    else
        limits!(ax, axesLimits[1], axesLimits[2], axesLimits[3], axesLimits[4], axesLimits[5], axesLimits[6])
    end

    # if (hideAll)
    #     hidedecorations!(ax)
    #     hidespines!(ax)
    # end

    if (showBox)
        hidedecorations!(ax, grid = false);
    else
        hidedecorations!(ax);
        hidespines!(ax);
    end


    fig
end

function ActiveFilaments.plotReachabilityCloudRGBSlice(sols, activationsGamma::Matrix{Float64}, gammaBounds, angle, sliceThickness, axesLimits; 
        gravity = false, full_solution = false, flipped = false, 
        resolution = (3840, 2160), 
        markersize = 2,
        transparency = false,
        opacity = 1.0,
        kwargs...)
    if full_solution
        if gravity
            points = [sol_i[end][1:3] for sol_i in sols]; # Assumign that second element of each sol is the support moment vector
        else
            points = sols.u; # ????
        end
    else
        if gravity
            points = [sol_i[1] for sol_i in sols]; # Assumign that second element of each sol is the support moment vector
        else
            points = sols.u;
        end
    end

    tana = (angle == pi / 2) ?  1 : tan(angle)
    fy = (angle == pi / 2) ?  0.0 : 1.0;
    t = (angle == pi / 2) ? sliceThickness : abs(sliceThickness / cos(angle));

    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    cropMap = ((fy * y - tana * x .<= t / 2) .& (fy * y - tana * x .>= -t / 2));
    xCrop = x .* cropMap;
    yCrop = y .* cropMap;
    zCrop = z .* cropMap;
    
    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    activations = [abs.(activationsGamma[i, :]) for i in axes(activationsGamma, 1)];
    
    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}();
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([gammaBounds[j][1][k], gammaBounds[j][2][k]]))));
            push!(gammaBoundsPairs, bounds);
        end
    end

    fig = Figure(resolution = resolution)
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = angle - pi/2, elevation = 0; kwargs...)
    GLMakie.scatter!(xCrop, yCrop, zCrop, 
        markersize = markersize, 
        color = [
            RGBA(activations[i][1] / gammaBoundsPairs[1][2], 
            activations[i][2] / gammaBoundsPairs[2][2], 
            activations[i][3] / gammaBoundsPairs[3][2], 
            opacity) 
            for i in 1:length(activations[1])
                ],
            transparency = transparency,
            perspectiveness = 0.0)

    hidedecorations!(ax)
    hidespines!(ax)
    if (flipped)
        limits!(ax, axesLimits[1], axesLimits[2], axesLimits[3], axesLimits[4], -axesLimits[6], -axesLimits[5])
    else
        limits!(ax, axesLimits[1], axesLimits[2], axesLimits[3], axesLimits[4], axesLimits[5], axesLimits[6])
    end

    fig
end

function ActiveFilaments.plotReachabilityCloudSlice(sols, activationsGamma::Matrix{Float64}, activationIndex::Integer, gammaBounds, angle, sliceThickness, axesLimits; 
        gravity = false, full_solution = false, flipped = false, 
        resolution = (3840, 2160), 
        markersize = 2,
        transparency = false,
        opacity = 1.0,
        colormap = :thermal,
        kwargs...)
    if full_solution
        if gravity
            points = [sol_i[end][1:3] for sol_i in sols]; # Assumign that second element of each sol is the support moment vector
        else
            points = sols.u; # ????
        end
    else
        if gravity
            points = [sol_i[1] for sol_i in sols]; # Assumign that second element of each sol is the support moment vector
        else
            points = sols.u;
        end
    end

    tana = (angle == pi / 2) ?  1 : tan(angle)
    fy = (angle == pi / 2) ?  0.0 : 1.0;
    t = (angle == pi / 2) ? sliceThickness : abs(sliceThickness / cos(angle));

    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    cropMap = ((fy * y - tana * x .<= t / 2) .& (fy * y - tana * x .>= -t / 2));
    xCrop = x .* cropMap;
    yCrop = y .* cropMap;
    zCrop = z .* cropMap;
    
    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    activations = [abs.(activationsGamma[:, i]) for i in axes(activationsGamma, 2)];

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}();
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([gammaBounds[j][1][k], gammaBounds[j][2][k]]))));
            push!(gammaBoundsPairs, bounds);
        end
    end

    fig = Figure(resolution = resolution)
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = angle - pi/2, elevation = 0; kwargs...)
    GLMakie.scatter!(xCrop, yCrop, zCrop, 
        markersize = markersize, 
        color = activations[activationIndex], 
        colormap = (colormap, opacity),
        transparency = transparency,
        perspectiveness = 0.0)

    hidedecorations!(ax)
    hidespines!(ax)
    if (flipped)
        limits!(ax, axesLimits[1], axesLimits[2], axesLimits[3], axesLimits[4], -axesLimits[6], -axesLimits[5])
    else
        limits!(ax, axesLimits[1], axesLimits[2], axesLimits[3], axesLimits[4], axesLimits[5], axesLimits[6])
    end

    fig
end


function computeTube(sol, R::AbstractFloat, Z, Θ; flipped = false)
    out = [sol(u)[1:3] + R * AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

function computeTube(sol, R::AbstractInterpolation, Z, Θ; flipped = false)
    out = [sol(u)[1:3] + R(u) * AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end


function computeTubeCap(sol, R, Z, Θ; flipped = false)
    out = [sol(Z)[1:3] + u * AngleAxis(v, sol(Z)[10], sol(Z)[11], sol(Z)[12]) * [sol(Z)[4], sol(Z)[5], sol(Z)[6]] for u in R, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

function computeFiberSurface(sol, R::AbstractFloat, R2::AbstractFloat, α2, Z, Θ; flipped = false)
    out = [sol(u)[1:3] + R * AngleAxis(v + u * tan(α2) / R2, sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

function computeFiberSurface(sol, R::AbstractInterpolation, Θ2T, Z, Θ; flipped = false)
    out = [sol(u)[1:3] + R(u) * AngleAxis(v + Θ2T(u), sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end


function computeFiberSide(sol, α2, Θ::AbstractFloat, Z, R; flipped = false)
    R2 = R[end];
    out = [sol(u)[1:3] + v * AngleAxis(Θ + u * tan(α2) / R2, sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in R]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

function computeFiberSide(sol, Θ, Θ2T, Z, R; flipped = false)
    out = Matrix{SVector{3, Float64}}(undef, length(Z), length(R(Z[1])));
    for i in eachindex(Z)
        u = Z[i];
        out[i, :] = [sol(u)[1:3] + v * AngleAxis(Θ + Θ2T(u), sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]] for v in R(u)]
        # for j in eachindex(R(Z[1]))
        #     u = Z[i];
        #     v = R(u);
        #     println(u);
        #     println(v);
        #     out[i, j] = sol(u)[1:3] + v * AngleAxis(Θ + Θ2T(u), sol(u)[10], sol(u)[11], sol(u)[12]) * [sol(u)[4], sol(u)[5], sol(u)[6]];
        # end
    end
    
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

function computeFiberCap(sol, α2::AbstractFloat, Z, Θ, R; flipped = false)
    R2 = R[end];
    out = [sol(Z)[1:3] + v * AngleAxis(u + Z * tan(α2) / R2, sol(Z)[10], sol(Z)[11], sol(Z)[12]) * [sol(Z)[4], sol(Z)[5], sol(Z)[6]] for u in Θ, v in R]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

function computeFiberCap(sol, Θ2T::Function, Z, Θ, R; flipped = false)
    out = [sol(Z)[1:3] + v * AngleAxis(u + Θ2T(Z), sol(Z)[10], sol(Z)[11], sol(Z)[12]) * [sol(Z)[4], sol(Z)[5], sol(Z)[6]] for u in Θ, v in R]
    x = getindex.(out, 1) * (flipped ? -1 : 1);
    y = getindex.(out, 2);
    z = getindex.(out, 3) * (flipped ? -1 : 1);

    return [x, y, z];
end

plotConfigurationTube!(filament::AFilament, sol; R_outer = filament.R0, flipped = false, n = 100, color = :black, opacity = 1.0) = 
    plotConfigurationTube!(Val(filament.tapered), filament, sol; R_outer = R_outer, flipped = flipped, n = n, color = color, opacity = opacity)

function plotConfigurationTube!(::Val{true}, filament, sol; R_outer = filament.R0, flipped = false, n = 100, color = :black, opacity = 1.0)
    Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n)

    if opacity < 1.0
        transparency = true;
    else
        transparency = false;
    end

    x, y, z = computeTube(sol, R_outer, Z, Θ; flipped = flipped)
    GLMakie.surface!(x, y, z, color = fill((color, opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false, transparency = transparency)

    R_tube_0 = range(0.0, R_outer(0.0), n);
    x, y, z = computeTubeCap(sol, R_tube_0, Z[1], Θ; flipped = flipped)
    GLMakie.surface!(x, y, z, color = fill((color, opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false, transparency = transparency)

    R_tube_end = range(0.0, R_outer(filament.L), n);
    x, y, z = computeTubeCap(sol, R_tube_end, Z[end], Θ; flipped = flipped)
    GLMakie.surface!(x, y, z, color = fill((color, opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false, transparency = transparency)
end

function plotConfigurationTube!(::Val{false}, filament, sol; R_outer = filament.R0, flipped = false, n = 100, color = :black, opacity = 1.0)
    Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n)

    if opacity < 1.0
        transparency = true;
    else
        transparency = false;
    end

    x, y, z = computeTube(sol, R_outer, Z, Θ; flipped = flipped)
    GLMakie.surface!(x, y, z, color = fill((color, opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false, transparency = transparency)

    R_tube = range(0.0, R_outer, n);
    x, y, z = computeTubeCap(sol, R_tube, Z[1], Θ; flipped = flipped)
    GLMakie.surface!(x, y, z, color = fill((color, opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false, transparency = transparency)

    x, y, z = computeTubeCap(sol, R_tube, Z[end], Θ; flipped = flipped)
    GLMakie.surface!(x, y, z, color = fill((color, opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false, transparency = transparency)
end

function plotConfigurationTubesSelfWeight!(filament, activation_structure, activations; 
    m0 = [0.0, 0.0, 0.0], uInit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, m0[1], m0[2], m0[3]], g = -9.8, Ng = 4, 
    kwargs...)

    g_range = range(start = 0.0, stop = g, length = Ng);

    for gamma in activations
        activation_conf = [copy(activation) for activation in activation_structure];
        for i in eachindex(gamma)
            activation_conf[i].γ = gamma[i];
        end

        activationFourier_conf = [piecewiseGammaToFourier(activation) for activation in activation_conf];
        
        sol_conf = selfWeightSolve(filament, activationFourier_conf, m0, uInit, g_range);

        plotConfigurationTube!(filament, sol_conf; kwargs...);
    end
end

function ActiveFilaments.plotFilamentCollapsedRings!(filament::AFilament{1}, activation_structure, sol; n = 100, flipped = false, opacity = 1.0, standard_core = true, perturb_deviation = 0.0, colors = [:yellow, :orange, :blue])
    #Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n) # Flip sign to flip normals

    Z = filament.Z;
    R_outer = standard_core ? cubic_spline_interpolation(Z, filament.rings[1].geometry.R2) : cubic_spline_interpolation(Z, filament.rings[1].geometry.R1);
    
    plotConfigurationTube!(filament, sol, R_outer = R_outer, flipped = flipped, n = n, color = :black, opacity = opacity);
    
    Z = range(0.0 - filament.L * perturb_deviation / 4.0, filament.L * (1 + perturb_deviation / 4), n);
    ZR = range(0.0 - filament.L * perturb_deviation / 4.0, filament.L * (1 + perturb_deviation / 4), length(filament.Z));
    for j in eachindex(filament.rings)
        θ0 = activation_structure[j].θ0;
        σ = activation_structure[j].σ;
        α2_r = filament.rings[j].fiberArchitecture.α2;
        R1_r = cubic_spline_interpolation(ZR, filament.rings[j].geometry.R1);
        R2_r = cubic_spline_interpolation(ZR, filament.rings[j].geometry.R2);
        R1_rp = cubic_spline_interpolation(ZR, filament.rings[j].geometry.R1 * (1 - perturb_deviation));
        R2_rp = cubic_spline_interpolation(ZR, filament.rings[j].geometry.R2 * (1 + perturb_deviation));
        Θ = range(θ0 - σ / 2, θ0 + σ / 2, n)
        Θ2T(Z) = -tan(α2_r) * log.(R2_r(Z) / filament.rings[j].geometry.R2[1]) / sin(filament.rings[j].geometry.phi2);
        R_range(Z) = range(R1_rp(Z), R2_rp(Z), n);

        x, y, z = computeFiberSurface(sol, R2_rp, Θ2T, Z, Θ; flipped = flipped)
        GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)

        x, y, z = computeFiberSurface(sol, R1_rp, Θ2T, Z, Θ; flipped = flipped)
        GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
        
        x, y, z = computeFiberSide(sol, Θ[1], Θ2T, Z, R_range; flipped = flipped)
        GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)

        x, y, z = computeFiberSide(sol, Θ[end], Θ2T, Z, R_range; flipped = flipped)
        GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
        
        x, y, z = computeFiberCap(sol, Θ2T, Z[1], Θ, R_range(Z[1]); flipped = flipped)
        GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
        
        x, y, z = computeFiberCap(sol, Θ2T, Z[end], Θ, R_range(Z[end]); flipped = flipped)
        GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
    end
end

function ActiveFilaments.plotFilamentCollapsedRings!(filament::AFilament{0}, activation_structure, sol; n = 100, flipped = false, opacity = 1.0, standard_core = true, perturb_deviation = 0.0, colors = [:yellow, :orange, :blue])
    Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n) # Flip sign to flip normals
    # set_theme!(background = false)
    # set_theme!()

    # R_tube = standard_core ? range(0.0, filament.rings[1].geometry.R2, n) : range(0.0, filament.rings[1].geometry.R1, n)
    R_outer = standard_core ? filament.rings[1].geometry.R2 : filament.rings[1].geometry.R1;
    
    plotConfigurationTube!(filament, sol, R_outer = R_outer, flipped = flipped, n = n, color = :black, opacity = opacity);
    
    α2 = filament.rings[1].fiberArchitecture.α2
    if (typeof(α2) == PiecewiseStructure)
        for i in eachindex(α2.piecewiseRanges)
            Zp = range(α2.piecewiseRanges[i][1] - (α2.piecewiseRanges[i][2] - α2.piecewiseRanges[i][1]) * perturb_deviation / 4,
                       α2.piecewiseRanges[i][2] + (α2.piecewiseRanges[i][2] - α2.piecewiseRanges[i][1]) * perturb_deviation / 4, n);
            for j in eachindex(filament.rings)
                θ0 = activation_structure[j].θ0;
                σ = activation_structure[j].σ;
                α2_r = filament.rings[j].fiberArchitecture.α2.expressions[i];
                R1_r = filament.rings[j].geometry.R1 * (1 - perturb_deviation);
                R2_r = filament.rings[j].geometry.R2 * (1 + perturb_deviation);
                Θ = range(θ0.expressions[i] - σ / 2, θ0.expressions[i] + σ / 2, n)
                R_range = range(R1_r, R2_r, n);
        
                x, y, z = computeFiberSurface(sol, R2_r, R2_r, α2_r, Zp, Θ; flipped = flipped)
                GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)

                x, y, z = computeFiberSurface(sol, R1_r, R2_r, α2_r, Zp, Θ; flipped = flipped)
                GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
                
                x, y, z = computeFiberSide(sol, α2_r, Θ[1], Z, R_range; flipped = flipped)
                GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)

                x, y, z = computeFiberSide(sol, α2_r, Θ[end], Z, R_range; flipped = flipped)
                GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
        
                x, y, z = computeFiberCap(sol, α2_r, Zp[1], Θ, R_range; flipped = flipped)
                GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
                
                x, y, z = computeFiberCap(sol, α2_r, Zp[end], Θ, R_range; flipped = flipped)
                GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
            end

            if standard_core
                # R2 = filament.rings[1].geometry.R2;
                # Θ_all = [];
                # θ0 = activation_structure[1].θ0;
                # σ = activation_structure[1].σ;
                # if (θ0.expressions[i] - σ / 2 > 0.0)
                #     append!(Θ_all, [range(0.0, θ0.expressions[i] - σ / 2, n)])
                # end
                # for j in 1:size(filament.rings) - 1
                #     θ0_1 = activation_structure[j].θ0;
                #     σ_1 = activation_structure[j].σ;
                #     θ0_2 = activation_structure[j + 1].θ0;
                #     σ_2 = activation_structure[j + 1].σ;
                #     append!(Θ_all, [range(θ0_1.expressions[i] + σ_1 / 2, θ0_2.expressions[i] - σ_2 / 2, n)]);
                # end
                # θ0 = activation_structure[end].θ0;
                # σ = activation_structure[end].σ;
                # append!(Θ_all, range(θ0.expressions[i] + σ / 2, 2 * pi, n))
                
                # x, y, z = computeTube(sol, R2, Zp, Θ; flipped = flipped)
            end
        end
    else
        Z = range(0.0 - filament.L * perturb_deviation / 4.0, filament.L * (1 + perturb_deviation / 4), n);
        for j in eachindex(filament.rings)
            θ0 = activation_structure[j].θ0;
            σ = activation_structure[j].σ;
            α2_r = filament.rings[j].fiberArchitecture.α2;
            R1_r = filament.rings[j].geometry.R1 * (1 - perturb_deviation);
            R2_r = filament.rings[j].geometry.R2 * (1 + perturb_deviation);
            Θ = range(θ0 - σ / 2, θ0 + σ / 2, n)
            R_range = range(R1_r, R2_r, n);
    
            x, y, z = computeFiberSurface(sol, R2_r, R2_r, α2_r, Z, Θ; flipped = flipped)
            GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)

            x, y, z = computeFiberSurface(sol, R1_r, R2_r, α2_r, Z, Θ; flipped = flipped)
            GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
            
            x, y, z = computeFiberSide(sol, α2_r, Θ[1], Z, R_range; flipped = flipped)
            GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)

            x, y, z = computeFiberSide(sol, α2_r, Θ[end], Z, R_range; flipped = flipped)
            GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
            
            x, y, z = computeFiberCap(sol, α2_r, Z[1], Θ, R_range; flipped = flipped)
            GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
            
            x, y, z = computeFiberCap(sol, α2_r, Z[end], Θ, R_range; flipped = flipped)
            GLMakie.surface!(x, y, z, color = fill((colors[j], opacity), n, n), shading = true, ssao = true, invert_normals = true, background = false)
        end
    end
end
end