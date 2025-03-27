#####################################################
### The main plotting module extension.           ###
### All visualization functions are defined here. ###
#####################################################
module PlottingExt

using ActiveFilaments
using ColorSchemes
using CairoMakie
using GLMakie
using Plots
using Interpolations
using Rotations
using StaticArrays
using DocStringExtensions

"""
    $(TYPEDSIGNATURES)

Plots the reachability cloud with RGB color-coding, given the
configuration solutions `sols`.
"""
function ActiveFilaments.plotReachabilityCloudRGB(sols, activationsGamma::Matrix{Float64},
        gammaBounds, axesLimits;
        gravity = false, full_solution = false, flipped = false, showBox = true,
        hideAll = false,
        azimuth = 1.275 * pi, elevation = pi / 8,
        resolution = (3840, 2160),
        size = nothing,
        markersize = 2,
        perspectiveness = 0.0,
        transparency = false,
        opacity = 1.0,
        for_rotation = false,
        size_to_ax = false,
        kwargs...)
    if full_solution
        if gravity
            # Assuming that the second element of each sol is the support moment vector
            points = [sol_i[end][1:3] for sol_i in sols]
        else
            # Case currently unused
        end
    else
        if gravity
            # Assuming that the second element of each sol is the support moment vector
            points = [sol_i[1] for sol_i in sols]
        else
            points = sols.u
        end
    end

    x = getindex.(points, 1) * (flipped ? -1 : 1)
    y = getindex.(points, 2)
    z = getindex.(points, 3) * (flipped ? -1 : 1)

    activations = [abs.(activationsGamma[i, :]) for i in axes(activationsGamma, 1)]

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}()
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([
                gammaBounds[j][1][k], gammaBounds[j][2][k]]))))
            push!(gammaBoundsPairs, bounds)
        end
    end

    GLMakie.activate!()
    if isnothing(size)
        fig = Figure(resolution = resolution)
    else
        fig = Figure(size = size)
    end
    if size_to_ax
        ax = Axis3(fig[1, 1], aspect = :data, azimuth = azimuth, elevation = elevation,
            perspectiveness = perspectiveness,
            width = size[1], height = size[2]; kwargs...)
    else
        ax = Axis3(
            fig[1, 1],
            aspect = :data,
            azimuth = azimuth,
            elevation = elevation,
            perspectiveness = perspectiveness;
            kwargs...
        )
    end
    if for_rotation
        ax.viewmode = :fit
    end
    GLMakie.scatter!(x, y, z, markersize = markersize,
        color = [RGBA(activations[i][1] / gammaBoundsPairs[1][2],
                     activations[i][2] / gammaBoundsPairs[2][2],
                     activations[i][3] / gammaBoundsPairs[3][2],
                     opacity)
                 for i in eachindex(activations)],
        transparency = transparency)

    if (flipped)
        limits!(
            ax,
            -axesLimits[2],
            -axesLimits[1],
            axesLimits[3],
            axesLimits[4],
            -axesLimits[6],
            -axesLimits[5]
        )
    else
        limits!(
            ax,
            axesLimits[1],
            axesLimits[2],
            axesLimits[3],
            axesLimits[4],
            axesLimits[5],
            axesLimits[6]
        )
    end

    if (showBox)
        hidedecorations!(ax, grid = false)
    else
        hidedecorations!(ax)
        hidespines!(ax)
    end

    fig
end

"""
    $(TYPEDSIGNATURES)

Plots the reachability cloud, given the
configuration solutions `sols`. The color-coding
is defined by `colormap`.
"""
function ActiveFilaments.plotReachabilityCloud(sols, activationsGamma::Matrix{Float64},
        activationIndex::Integer, gammaBounds, axesLimits;
        gravity = false, full_solution = false, flipped = false, showBox = true,
        hideAll = false,
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
            # Assuming that second element of each sol is the support moment vector
            points = [sol_i[end][1:3] for sol_i in sols]
        else
            # Case currently unused
        end
    else
        if gravity
            # Assuming that second element of each sol is the support moment vector
            points = [sol_i[1] for sol_i in sols]
        else
            points = sols.u
        end
    end

    x = getindex.(points, 1) * (flipped ? -1 : 1)
    y = getindex.(points, 2)
    z = getindex.(points, 3) * (flipped ? -1 : 1)

    activations = [abs.(activationsGamma[:, i]) for i in axes(activationsGamma, 2)]

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}()
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([
                gammaBounds[j][1][k], gammaBounds[j][2][k]]))))
            push!(gammaBoundsPairs, bounds)
        end
    end

    if !isnothing(maxAbsGamma)
        if (maxAbsGammaIndexed)
            activationMap = activations[activationIndex] .<= maxAbsGamma
        else
            activationMap = (activations[1] .<= maxAbsGamma)
            for i in 2:size(activationsGamma, 2)
                activationMap = activationMap .& (activations[i] .<= maxAbsGamma)
            end
        end

        x = x[activationMap]
        y = y[activationMap]
        z = z[activationMap]

        for i in axes(activationsGamma, 2)
            activations[i] = activations[i][activationMap]
        end
    end

    GLMakie.activate!()
    fig = Figure(resolution = resolution)
    ax = Axis3(
        fig[1, 1],
        aspect = :data,
        azimuth = azimuth,
        elevation = elevation,
        perspectiveness = perspectiveness;
        kwargs...
    )
    GLMakie.scatter!(
        x,
        y,
        z,
        markersize = markersize,
        color = activations[activationIndex],
        colormap = (colormap, opacity),
        colorrange = gammaBoundsPairs[activationIndex],
        transparency = transparency
    )

    if (flipped)
        limits!(
            ax,
            -axesLimits[2],
            -axesLimits[1],
            axesLimits[3],
            axesLimits[4],
            -axesLimits[6],
            -axesLimits[5]
        )
    else
        limits!(
            ax,
            axesLimits[1],
            axesLimits[2],
            axesLimits[3],
            axesLimits[4],
            axesLimits[5],
            axesLimits[6]
        )
    end

    if (showBox)
        hidedecorations!(ax, grid = false)
    else
        hidedecorations!(ax)
        hidespines!(ax)
    end

    fig
end

"""
    $(TYPEDSIGNATURES)

Plots the reachability cloud slice with RGB color-coding, given the
configuration solutions `sols`, a slice angle `angle`, and a slice thickness
`sliceThickness`.
"""
function ActiveFilaments.plotReachabilityCloudRGBSlice(sols,
        activationsGamma::Matrix{Float64}, gammaBounds, angle, sliceThickness, axesLimits;
        gravity = false, full_solution = false, flipped = false,
        resolution = (3840, 2160),
        markersize = 2,
        transparency = false,
        opacity = 1.0,
        kwargs...)
    if full_solution
        if gravity
            # Assuming that second element of each sol is the support moment vector
            points = [sol_i[end][1:3] for sol_i in sols]
        else
            # Case currently unused
        end
    else
        if gravity
            # Assuming that second element of each sol is the support moment vector
            points = [sol_i[1] for sol_i in sols]
        else
            points = sols.u
        end
    end

    tana = (angle == pi / 2) ? 1 : tan(angle)
    fy = (angle == pi / 2) ? 0.0 : 1.0
    t = (angle == pi / 2) ? sliceThickness : abs(sliceThickness / cos(angle))

    x = getindex.(points, 1) * (flipped ? -1 : 1)
    y = getindex.(points, 2)
    z = getindex.(points, 3) * (flipped ? -1 : 1)

    cropMap = ((fy * y - tana * x .<= t / 2) .& (fy * y - tana * x .>= -t / 2))

    xCrop = x[cropMap]
    yCrop = y[cropMap]
    zCrop = z[cropMap]

    activations = [abs.(activationsGamma[i, :]) for i in axes(activationsGamma, 1)]

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}()
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([
                gammaBounds[j][1][k], gammaBounds[j][2][k]]))))
            push!(gammaBoundsPairs, bounds)
        end
    end

    fig = Figure(resolution = resolution)
    ax = Axis3(
        fig[1, 1], aspect = :data, azimuth = angle - pi / 2, elevation = 0; kwargs...)
    GLMakie.scatter!(xCrop, yCrop, zCrop,
        markersize = markersize,
        color = [RGBA(activations[i][1] / gammaBoundsPairs[1][2],
                     activations[i][2] / gammaBoundsPairs[2][2],
                     activations[i][3] / gammaBoundsPairs[3][2],
                     opacity)
                 for i in 1:length(activations[1])],
        transparency = transparency,
        perspectiveness = 0.0)

    hidedecorations!(ax)
    hidespines!(ax)
    if (flipped)
        limits!(
            ax,
            axesLimits[1],
            axesLimits[2],
            axesLimits[3],
            axesLimits[4],
            -axesLimits[6],
            -axesLimits[5]
        )
    else
        limits!(
            ax,
            axesLimits[1],
            axesLimits[2],
            axesLimits[3],
            axesLimits[4],
            axesLimits[5],
            axesLimits[6]
        )
    end

    fig
end

"""
    $(TYPEDSIGNATURES)

Plots the reachability cloud slice, given the
configuration solutions `sols`, a slice angle `angle`, and a slice thickness
`sliceThickness`. The color-coding is defined by `colormap`.
"""
function ActiveFilaments.plotReachabilityCloudSlice(
        sols, activationsGamma::Matrix{Float64},
        activationIndex::Integer, gammaBounds, angle, sliceThickness, axesLimits;
        gravity = false, full_solution = false, flipped = false,
        resolution = (3840, 2160),
        markersize = 2,
        transparency = false,
        opacity = 1.0,
        colormap = :thermal,
        kwargs...)
    if full_solution
        if gravity
            # Assuming that second element of each sol is the support moment vector
            points = [sol_i[end][1:3] for sol_i in sols]
        else
            # Case currently unused
        end
    else
        if gravity
            # Assuming that second element of each sol is the support moment vector
            points = [sol_i[1] for sol_i in sols]
        else
            points = sols.u
        end
    end

    tana = (angle == pi / 2) ? 1 : tan(angle)
    fy = (angle == pi / 2) ? 0.0 : 1.0
    t = (angle == pi / 2) ? sliceThickness : abs(sliceThickness / cos(angle))

    x = getindex.(points, 1) * (flipped ? -1 : 1)
    y = getindex.(points, 2)
    z = getindex.(points, 3) * (flipped ? -1 : 1)

    cropMap = ((fy * y - tana * x .<= t / 2) .& (fy * y - tana * x .>= -t / 2))

    xCrop = x[cropMap]
    yCrop = y[cropMap]
    zCrop = z[cropMap]

    activations = [abs.(activationsGamma[:, i]) for i in axes(activationsGamma, 2)]

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}()
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([
                gammaBounds[j][1][k], gammaBounds[j][2][k]]))))
            push!(gammaBoundsPairs, bounds)
        end
    end

    fig = Figure(resolution = resolution)
    ax = Axis3(
        fig[1, 1], aspect = :data, azimuth = angle - pi / 2, elevation = 0; kwargs...)
    GLMakie.scatter!(xCrop, yCrop, zCrop,
        markersize = markersize,
        color = activations[activationIndex][cropMap],
        colormap = (colormap, opacity),
        transparency = transparency,
        perspectiveness = 0.0)

    hidedecorations!(ax)
    hidespines!(ax)
    if (flipped)
        limits!(
            ax,
            axesLimits[1],
            axesLimits[2],
            axesLimits[3],
            axesLimits[4],
            -axesLimits[6],
            -axesLimits[5]
        )
    else
        limits!(
            ax,
            axesLimits[1],
            axesLimits[2],
            axesLimits[3],
            axesLimits[4],
            axesLimits[5],
            axesLimits[6]
        )
    end

    fig
end

function computeTube(sol, R::AbstractFloat, Z, Θ; flipped = false)
    out = [sol(u)[1:3] +
           R * AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeTube(sol, R::AbstractInterpolation, Z, Θ; flipped = false)
    out = [sol(u)[1:3] +
           R(u) * AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeTubeCap(sol, R, Z, Θ; flipped = false)
    out = [sol(Z)[1:3] +
           u * AngleAxis(v, sol(Z)[10], sol(Z)[11], sol(Z)[12]) *
           [sol(Z)[4], sol(Z)[5], sol(Z)[6]] for u in R, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeFiberSurface(
        sol,
        R::AbstractFloat,
        R2::AbstractFloat,
        α2,
        Z,
        Θ;
        flipped = false
)
    out = [sol(u)[1:3] +
           R * AngleAxis(v + u * tan(α2) / R2, sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeFiberSurface(sol, R::AbstractInterpolation, Θ2T, Z, Θ; flipped = false)
    out = [sol(u)[1:3] +
           R(u) * AngleAxis(v + Θ2T(u), sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeFiberSide(sol, α2, Θ::AbstractFloat, Z, R; flipped = false)
    R2 = R[end]
    out = [sol(u)[1:3] +
           v * AngleAxis(Θ + u * tan(α2) / R2, sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in R]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeFiberSide(sol, Θ, Θ2T, Z, R; flipped = false)
    out = Matrix{SVector{3, Float64}}(undef, length(Z), length(R(Z[1])))
    for i in eachindex(Z)
        u = Z[i]
        out[i, :] = [sol(u)[1:3] +
                     v * AngleAxis(Θ + Θ2T(u), sol(u)[10], sol(u)[11], sol(u)[12]) *
                     [sol(u)[4], sol(u)[5], sol(u)[6]] for v in R(u)]
    end

    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeFiberCap(sol, α2::AbstractFloat, Z, Θ, R; flipped = false)
    R2 = R[end]
    out = [sol(Z)[1:3] +
           v * AngleAxis(u + Z * tan(α2) / R2, sol(Z)[10], sol(Z)[11], sol(Z)[12]) *
           [sol(Z)[4], sol(Z)[5], sol(Z)[6]] for u in Θ, v in R]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function computeFiberCap(sol, Θ2T::Function, Z, Θ, R; flipped = false)
    out = [sol(Z)[1:3] +
           v * AngleAxis(u + Θ2T(Z), sol(Z)[10], sol(Z)[11], sol(Z)[12]) *
           [sol(Z)[4], sol(Z)[5], sol(Z)[6]] for u in Θ, v in R]
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return [x, y, z]
end

function plotConfigurationTube!(
        filament::AFilament,
        sol;
        R_outer = filament.tapered ? 
                  cubic_spline_interpolation(filament.Z, filament.R0) : 
                  filament.R0,
        flipped = false,
        n = 100,
        color = :black,
        opacity = 1.0
)
    plotConfigurationTube!(
        Val(filament.tapered),
        filament,
        sol;
        R_outer = R_outer,
        flipped = flipped,
        n = n,
        color = color,
        opacity = opacity
    )
end

function plotConfigurationTube!(
        ::Val{true},
        filament,
        sol;
        R_outer = filament.R0,
        flipped = false,
        n = 100,
        color = :black,
        opacity = 1.0
)
    Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n)

    if opacity < 1.0
        transparency = true
    else
        transparency = false
    end

    x, y, z = computeTube(sol, R_outer, Z, Θ; flipped = flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        shading = Makie.automatic,
        ssao = true,
        invert_normals = true,
        background = false,
        transparency = transparency
    )

    R_tube_0 = range(0.0, R_outer(0.0), n)
    x, y, z = computeTubeCap(sol, R_tube_0, Z[1], Θ; flipped = flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        shading = Makie.automatic,
        ssao = true,
        invert_normals = true,
        background = false,
        transparency = transparency
    )

    R_tube_end = range(0.0, R_outer(filament.L), n)
    x, y, z = computeTubeCap(sol, R_tube_end, Z[end], Θ; flipped = flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        shading = Makie.automatic,
        ssao = true,
        invert_normals = true,
        background = false,
        transparency = transparency
    )
end

function plotConfigurationTube!(
        ::Val{false},
        filament,
        sol;
        R_outer = filament.R0,
        flipped = false,
        n = 100,
        color = :black,
        opacity = 1.0
)
    Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n)

    if opacity < 1.0
        transparency = true
    else
        transparency = false
    end

    x, y, z = computeTube(sol, R_outer, Z, Θ; flipped = flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        shading = Makie.automatic,
        ssao = true,
        invert_normals = true,
        background = false,
        transparency = transparency
    )

    R_tube = range(0.0, R_outer, n)
    x, y, z = computeTubeCap(sol, R_tube, Z[1], Θ; flipped = flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        shading = Makie.automatic,
        ssao = true,
        invert_normals = true,
        background = false,
        transparency = transparency
    )

    x, y, z = computeTubeCap(sol, R_tube, Z[end], Θ; flipped = flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        shading = Makie.automatic,
        ssao = true,
        invert_normals = true,
        background = false,
        transparency = transparency
    )
end

function plotConfigurationTubesSelfWeight!(filament, activation_structure, activations;
        bcs = [
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
            ],
        m0 = [0.0, 0.0, 0.0],
        u0 = [bcs..., m0...],
        g = -9.8, 
        Ng = 2,
        solver = 1,
        kwargs...)
    g_range = range(start = 0.0, stop = g, length = Ng)

    for gamma in activations
        activation_conf = [copy(activation) for activation in activation_structure]
        for i in eachindex(gamma)
            activation_conf[i].γ = gamma[i]
        end

        activationFourier_conf = [piecewiseGammaToFourier(activation)
                                  for activation in activation_conf]

        sol_conf = selfWeightSolve(filament, activationFourier_conf, u0, bcs, g_range; solver = solver)

        plotConfigurationTube!(filament, sol_conf; kwargs...)
    end
end

"""
    $(TYPEDSIGNATURES)

Plots the `filament` with an activation format `activation_structure`,
and a configuration solution `sol` for the case of multiple overlapping
rings.
"""
function ActiveFilaments.plotFilamentCollapsedRings!(
        filament::AFilament{1},
        activation_structure,
        sol;
        n = 100,
        flipped = false,
        opacity_core = 1.0,
        opacity_fibers = 1.0,
        standard_core = true,
        perturb_deviation = 0.0,
        colors = [:yellow, :orange, :blue]
)
    Θ = range(0, 2 * pi, n)

    Z = filament.Z
    R_outer = standard_core ? cubic_spline_interpolation(Z, filament.rings[1].geometry.R2) :
              cubic_spline_interpolation(Z, filament.rings[1].geometry.R1)

    plotConfigurationTube!(
        filament,
        sol,
        R_outer = R_outer,
        flipped = flipped,
        n = n,
        color = :black,
        opacity = opacity_core
    )

    Z = range(
        0.0 - filament.L * perturb_deviation / 4.0,
        filament.L * (1 + perturb_deviation / 4),
        n
    )
    ZR = range(
        0.0 - filament.L * perturb_deviation / 4.0,
        filament.L * (1 + perturb_deviation / 4),
        length(filament.Z)
    )
    for j in eachindex(filament.rings)
        θ0 = activation_structure[j].θ0
        σ = activation_structure[j].σ
        α2_r = filament.rings[j].fiberArchitecture.α2
        R1_r = cubic_spline_interpolation(ZR, filament.rings[j].geometry.R1)
        R2_r = cubic_spline_interpolation(ZR, filament.rings[j].geometry.R2)
        R1_rp = cubic_spline_interpolation(
            ZR,
            filament.rings[j].geometry.R1 * (1 - perturb_deviation)
        )
        R2_rp = cubic_spline_interpolation(
            ZR,
            filament.rings[j].geometry.R2 * (1 + perturb_deviation)
        )
        Θ = range(θ0 - σ / 2, θ0 + σ / 2, n)
        function Θ2T(Z)
            -tan(α2_r) * log.(R2_r(Z) / filament.rings[j].geometry.R2[1]) /
            sin(filament.rings[j].geometry.phi2)
        end
        R_range(Z) = range(R1_rp(Z), R2_rp(Z), n)

        x, y, z = computeFiberSurface(sol, R2_rp, Θ2T, Z, Θ; flipped = flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((colors[j], opacity_fibers), n, n),
            shading = Makie.Automatic(),
            ssao = true,
            invert_normals = true,
            background = false
        )

        x, y, z = computeFiberSurface(sol, R1_rp, Θ2T, Z, Θ; flipped = flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((colors[j], opacity_fibers), n, n),
            shading = Makie.Automatic(),
            ssao = true,
            invert_normals = true,
            background = false
        )

        x, y, z = computeFiberSide(sol, Θ[1], Θ2T, Z, R_range; flipped = flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((colors[j], opacity_fibers), n, n),
            shading = Makie.Automatic(),
            ssao = true,
            invert_normals = true,
            background = false
        )

        x, y, z = computeFiberSide(sol, Θ[end], Θ2T, Z, R_range; flipped = flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((colors[j], opacity_fibers), n, n),
            shading = Makie.Automatic(),
            ssao = true,
            invert_normals = true,
            background = false
        )

        x, y, z = computeFiberCap(sol, Θ2T, Z[1], Θ, R_range(Z[1]); flipped = flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((colors[j], opacity_fibers), n, n),
            shading = Makie.Automatic(),
            ssao = true,
            invert_normals = true,
            background = false
        )

        x, y, z = computeFiberCap(sol, Θ2T, Z[end], Θ, R_range(Z[end]); flipped = flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((colors[j], opacity_fibers), n, n),
            shading = Makie.Automatic(),
            ssao = true,
            invert_normals = true,
            background = false
        )
    end
end

function ActiveFilaments.plotFilamentCollapsedRings!(filament::AFilament{0},
        activation_structure, sol;
        n = 100, flipped = false, opacity_core = 1.0, opacity_fibers = 1.0,
        standard_core = true, perturb_deviation = 0.0, colors = [:yellow, :orange, :blue])
    Z = range(0, filament.L, n)
    Θ = range(0, 2 * pi, n) # Flip sign to flip normals

    R_outer = standard_core ? filament.rings[1].geometry.R2 : filament.rings[1].geometry.R1

    plotConfigurationTube!(
        filament,
        sol,
        R_outer = R_outer,
        flipped = flipped,
        n = n,
        color = :black,
        opacity = opacity_core
    )

    α2 = filament.rings[1].fiberArchitecture.α2
    if (typeof(α2) == PiecewiseStructure)
        for i in eachindex(α2.piecewiseRanges)
            Zp = range(
                α2.piecewiseRanges[i][1] -
                (α2.piecewiseRanges[i][2] - α2.piecewiseRanges[i][1]) * perturb_deviation /
                4,
                α2.piecewiseRanges[i][2] +
                (α2.piecewiseRanges[i][2] - α2.piecewiseRanges[i][1]) * perturb_deviation /
                4, n)
            for j in eachindex(filament.rings)
                θ0 = activation_structure[j].θ0
                σ = activation_structure[j].σ
                α2_r = filament.rings[j].fiberArchitecture.α2.expressions[i]
                R1_r = filament.rings[j].geometry.R1 * (1 - perturb_deviation)
                R2_r = filament.rings[j].geometry.R2 * (1 + perturb_deviation)
                Θ = range(θ0.expressions[i] - σ / 2, θ0.expressions[i] + σ / 2, n)
                R_range = range(R1_r, R2_r, n)

                x, y, z = computeFiberSurface(
                    sol, R2_r, R2_r, α2_r, Zp, Θ; flipped = flipped)
                GLMakie.surface!(
                    x,
                    y,
                    z,
                    color = fill((colors[j], opacity_fibers), n, n),
                    shading = Makie.Automatic(),
                    ssao = true,
                    invert_normals = true,
                    background = false
                )

                x, y, z = computeFiberSurface(
                    sol, R1_r, R2_r, α2_r, Zp, Θ; flipped = flipped)
                GLMakie.surface!(
                    x,
                    y,
                    z,
                    color = fill((colors[j], opacity_fibers), n, n),
                    shading = Makie.Automatic(),
                    ssao = true,
                    invert_normals = true,
                    background = false
                )

                x, y, z = computeFiberSide(sol, α2_r, Θ[1], Z, R_range; flipped = flipped)
                GLMakie.surface!(
                    x,
                    y,
                    z,
                    color = fill((colors[j], opacity_fibers), n, n),
                    shading = Makie.Automatic(),
                    ssao = true,
                    invert_normals = true,
                    background = false
                )

                x, y, z = computeFiberSide(sol, α2_r, Θ[end], Z, R_range; flipped = flipped)
                GLMakie.surface!(
                    x,
                    y,
                    z,
                    color = fill((colors[j], opacity_fibers), n, n),
                    shading = Makie.Automatic(),
                    ssao = true,
                    invert_normals = true,
                    background = false
                )

                x, y, z = computeFiberCap(sol, α2_r, Zp[1], Θ, R_range; flipped = flipped)
                GLMakie.surface!(
                    x,
                    y,
                    z,
                    color = fill((colors[j], opacity_fibers), n, n),
                    shading = Makie.Automatic(),
                    ssao = true,
                    invert_normals = true,
                    background = false
                )

                x, y, z = computeFiberCap(sol, α2_r, Zp[end], Θ, R_range; flipped = flipped)
                GLMakie.surface!(
                    x,
                    y,
                    z,
                    color = fill((colors[j], opacity_fibers), n, n),
                    shading = Makie.Automatic(),
                    ssao = true,
                    invert_normals = true,
                    background = false
                )
            end
        end
    else
        Z = range(
            0.0 - filament.L * perturb_deviation / 4.0,
            filament.L * (1 + perturb_deviation / 4),
            n
        )
        for j in eachindex(filament.rings)
            θ0 = activation_structure[j].θ0
            σ = activation_structure[j].σ
            α2_r = filament.rings[j].fiberArchitecture.α2
            R1_r = filament.rings[j].geometry.R1 * (1 - perturb_deviation)
            R2_r = filament.rings[j].geometry.R2 * (1 + perturb_deviation)
            Θ = range(θ0 - σ / 2, θ0 + σ / 2, n)
            R_range = range(R1_r, R2_r, n)

            x, y, z = computeFiberSurface(sol, R2_r, R2_r, α2_r, Z, Θ; flipped = flipped)
            GLMakie.surface!(
                x,
                y,
                z,
                color = fill((colors[j], opacity_fibers), n, n),
                shading = Makie.Automatic(),
                invert_normals = true,
                ssao = true,
                background = false
            )

            x, y, z = computeFiberSurface(sol, R1_r, R2_r, α2_r, Z, Θ; flipped = flipped)
            GLMakie.surface!(
                x,
                y,
                z,
                color = fill((colors[j], opacity_fibers), n, n),
                shading = Makie.Automatic(),
                invert_normals = true,
                ssao = true,
                background = false
            )

            x, y, z = computeFiberSide(sol, α2_r, Θ[1], Z, R_range; flipped = flipped)
            GLMakie.surface!(
                x,
                y,
                z,
                color = fill((colors[j], opacity_fibers), n, n),
                shading = Makie.Automatic(),
                invert_normals = true,
                ssao = true,
                background = false
            )

            x, y, z = computeFiberSide(sol, α2_r, Θ[end], Z, R_range; flipped = flipped)
            GLMakie.surface!(
                x,
                y,
                z,
                color = fill((colors[j], opacity_fibers), n, n),
                shading = Makie.Automatic(),
                invert_normals = true,
                ssao = true,
                background = false
            )

            x, y, z = computeFiberCap(sol, α2_r, Z[1], Θ, R_range; flipped = flipped)
            GLMakie.surface!(
                x,
                y,
                z,
                color = fill((colors[j], opacity_fibers), n, n),
                shading = Makie.Automatic(),
                invert_normals = true,
                ssao = true,
                background = false
            )

            x, y, z = computeFiberCap(sol, α2_r, Z[end], Θ, R_range; flipped = flipped)
            GLMakie.surface!(
                x,
                y,
                z,
                color = fill((colors[j], opacity_fibers), n, n),
                shading = Makie.Automatic(),
                invert_normals = true,
                ssao = true,
                background = false
            )
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Plots the `filament` with an activation format `activation_structure`,
for all sets of fibrillar activations provided in `activations` assuming
the self-weight loading scenario. This plotting function visualizes a simplified
tubular representation of the deformed `filament`.
"""
function ActiveFilaments.plotConfigurationTubesSelfWeight!(filament, activation_structure,
        activations;
        bcs = [
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
            ],
        m0 = [0.0, 0.0, 0.0],
        u0 = [bcs..., m0...],
        g = -9.8, 
        Ng = 2,
        solver = 1,
        kwargs...)
    g_range = range(start = 0.0, stop = g, length = Ng)

    for gamma in activations
        activation_conf = [copy(activation) for activation in activation_structure]
        for i in eachindex(gamma)
            activation_conf[i].γ = gamma[i]
        end

        activationFourier_conf::Vector{ActivationFourier} = [piecewiseGammaToFourier(activation)
                                                             for activation in activation_conf]

        sol_conf = selfWeightSolve(filament, activationFourier_conf, u0, bcs, g_range; solver = solver)

        plotConfigurationTube!(filament, sol_conf; kwargs...)
    end
end

"""
    $(TYPEDSIGNATURES)

Plots the `filament` with an activation format `activation_structure`,
for all sets of fibrillar activations provided in `activations` assuming
the self-weight loading scenario.
"""
function ActiveFilaments.plotConfigurationsSelfWeight!(
        filament,
        activation_structure,
        activations;
        bcs = [
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
            ],
        m0 = [0.0, 0.0, 0.0],
        u0 = [bcs..., m0...],
        g = -9.8,
        Ng = 2,
        solver = 1,
        kwargs...
)
    g_range = range(start = 0.0, stop = g, length = Ng)

    for gamma in activations
        activation_conf = [copy(activation) for activation in activation_structure]
        for i in eachindex(gamma)
            activation_conf[i].γ = gamma[i]
        end

        activationFourier_conf::Vector{ActivationFourier} = [piecewiseGammaToFourier(activation)
                                                             for activation in activation_conf]

        sol_conf = selfWeightSolve(filament, activationFourier_conf, u0, bcs, g_range; solver = solver)

        plotFilamentCollapsedRings!(filament, activation_structure, sol_conf; kwargs...)
    end
end

##########################################################
### Plotting functions for the Trunk package expansion ###
##########################################################
function extract_xyz(out, flipped)
    x = getindex.(out, 1) * (flipped ? -1 : 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3) * (flipped ? -1 : 1)

    return x, y, z
end

function plot_muscle!(
        sol,
        R1,
        R2,
        Θ1,
        Θ2,
        Z1,
        Z2,
        color,
        R_factor = 1.0;
        n = 40,
        flipped = false,
        opacity_fibers = 1.0,
        R1_surface = false,
        Z1_cap = false,
        transparency = false,
        override_transparent = false
)
    Θ = LinRange(Θ1, Θ2, n)
    Z = LinRange(Z1, Z2, n)

    # R2 surface
    out = [sol(u)[1:3] +
           R_factor * R2(u) * AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity_fibers), n, n),
        ssao = true,
        invert_normals = true,
        transparency = transparency
    )

    # R1 surface
    if R1_surface
        out = [sol(u)[1:3] +
               R_factor * R1(u) * AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) *
               [sol(u)[4], sol(u)[5], sol(u)[6]] for u in Z, v in Θ]
        x, y, z = extract_xyz(out, flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((color, opacity_fibers), n, n),
            ssao = true,
            invert_normals = false,
            transparency = transparency
        )
    end

    # Z1 cap
    if Z1_cap
        R = LinRange(R_factor * R1(Z1), R_factor * R2(Z1), n)
        out = [sol(Z1)[1:3] +
               v * AngleAxis(u, sol(Z1)[10], sol(Z1)[11], sol(Z1)[12]) *
               [sol(Z1)[4], sol(Z1)[5], sol(Z1)[6]] for u in Θ, v in R]
        x, y, z = extract_xyz(out, flipped)
        GLMakie.surface!(
            x,
            y,
            z,
            color = fill((color, override_transparent ? opacity_fibers : 1.0), n, n),
            ssao = true,
            invert_normals = false,
            transparency = override_transparent
        )
    end

    # Z2 cap
    R = LinRange(R_factor * R1(Z2), R_factor * R2(Z2), n)
    out = [sol(Z2)[1:3] +
           v * AngleAxis(u, sol(Z2)[10], sol(Z2)[11], sol(Z2)[12]) *
           [sol(Z2)[4], sol(Z2)[5], sol(Z2)[6]] for u in Θ, v in R]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, override_transparent ? opacity_fibers : 1.0), n, n),
        ssao = true,
        invert_normals = true,
        transparency = override_transparent
    )

    # Θ1 side
    out = Matrix{SVector{3, Float64}}(undef, n, n)
    for i in eachindex(Z)
        u = Z[i]
        out[i, :] = [sol(u)[1:3] +
                     v * AngleAxis(Θ1, sol(u)[10], sol(u)[11], sol(u)[12]) *
                     [sol(u)[4], sol(u)[5], sol(u)[6]]
                     for
                     v in LinRange(R_factor * R1(u), R_factor * R2(u), n)]
    end
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity_fibers), n, n),
        ssao = true,
        invert_normals = true,
        transparency = transparency
    )

    # Θ2 side
    out = Matrix{SVector{3, Float64}}(undef, n, n)
    for i in eachindex(Z)
        u = Z[i]
        out[i, :] = [sol(u)[1:3] +
                     v * AngleAxis(Θ2, sol(u)[10], sol(u)[11], sol(u)[12]) *
                     [sol(u)[4], sol(u)[5], sol(u)[6]]
                     for
                     v in LinRange(R_factor * R1(u), R_factor * R2(u), n)]
    end
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity_fibers), n, n),
        ssao = true,
        invert_normals = false,
        transparency = transparency
    )
end

function plot_outer_trunk!(
        trunk::Trunk{T, N},
        sol,
        R_factor = 1.0;
        n = 40,
        color = :gray,
        opacity = 0.3,
        flipped = false
) where {T, N}
    out = [sol(u)[1:3] +
           (R_factor * (trunk.R00 * 1.01 - u * tan(trunk.φ0))) *
           AngleAxis(v, sol(u)[10], sol(u)[11], sol(u)[12]) *
           [sol(u)[4], sol(u)[5], sol(u)[6]]
           for u in LinRange(0.0, trunk.L, n),
    v in LinRange(0.0, 2 * pi, n)]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = true,
        transparency = true
    )

    out = [sol(0.0)[1:3] - trunk.L / 1000.0 * sol(0.0)[10:12] +
           v * AngleAxis(u, sol(0.0)[10], sol(0.0)[11], sol(0.0)[12]) *
           [sol(0.0)[4], sol(0.0)[5], sol(0.0)[6]]
           for u in LinRange(0.0, 2 * pi, n),
    v in LinRange(0.0, R_factor * trunk.R00 * 1.01, n)]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = false,
        transparency = true
    )

    out = [sol(trunk.L)[1:3] + trunk.L / 1000.0 * sol(trunk.L)[10:12] +
           v * AngleAxis(u, sol(trunk.L)[10], sol(trunk.L)[11], sol(trunk.L)[12]) *
           [sol(trunk.L)[4], sol(trunk.L)[5], sol(trunk.L)[6]]
           for
           u in LinRange(0.0, 2 * pi, n),
    v in LinRange(0.0, R_factor * (trunk.R00 * 1.01 - trunk.L * tan(trunk.φ0)), n)]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = true,
        transparency = true
    )
end

"""
    $(TYPEDSIGNATURES)

Plots a standard visualization of the trunk `trunk_sim` 
for a configuration solution `sol`, and the auxiliary 
struct `a` of type `ActivatedTrunkQuantities`.
"""
function ActiveFilaments.plot_trunk!(
        trunk_sim::TrunkFast{T, N},
        sol,
        a;
        n = 40,
        colors = [:orange, :magenta, :cyan, :red, :blue],
        flipped = false,
        opacity_fibers = 1.0,
        opacity_trunk = 0.3,
        override_transparent = false
) where {T, N}
    trunk = trunk_sim.trunk
    R1 = trunk_sim.interpolations.R1
    R2 = trunk_sim.interpolations.R2

    R_factor = compute_R_factor_current(trunk_sim.trunk, sol)
    
    plot_outer_trunk!(
        trunk,
        sol,
        R_factor;
        n = n,
        flipped = flipped,
        opacity = opacity_trunk
    )
    for i in 1:T
        for j in 1:5
            opacity_f = override_transparent ? opacity_fibers :
                        (j < 5 ? opacity_fibers : 1.0)
            transp = override_transparent ? true : (j < 5)
            plot_muscle!(
                sol,
                R1[j](i),
                R2[j](i),
                trunk.Θ1R[i, j],
                trunk.Θ2R[i, j],
                trunk.Z1[i],
                trunk.Z2[i],
                colors[j],
                R_factor;
                n = n,
                R1_surface = (j >= 4),
                Z1_cap = (i == 1),
                flipped = flipped,
                opacity_fibers = opacity_f,
                transparency = transp,
                override_transparent = override_transparent
            )

            plot_muscle!(
                sol,
                R1[j](i),
                R2[j](i),
                trunk.Θ1L[i, j],
                trunk.Θ2L[i, j],
                trunk.Z1[i],
                trunk.Z2[i],
                colors[j],
                R_factor;
                n = n,
                R1_surface = (j >= 4),
                Z1_cap = (i == 1),
                flipped = flipped,
                opacity_fibers = opacity_f,
                transparency = transp,
                override_transparent = override_transparent
            )
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Plots an isolated set of muscles of the trunk `trunk_sim` 
for a configuration solution `sol`, the auxiliary struct
`a` of type `ActivatedTrunkQuantities`,
the indices `muscle_indices` of the muscles to be plotted, and
the `colors` of the muscle surfaces.
"""
function ActiveFilaments.plot_trunk_isolated!(
        trunk_sim::TrunkFast{T, N},
        sol,
        a,
        muscle_indices,
        colors;
        n = 40,
        flipped = false,
        opacity_fibers = 1.0,
        opacity_trunk = 0.3
) where {T, N}
    trunk = trunk_sim.trunk
    R1 = trunk_sim.interpolations.R1
    R2 = trunk_sim.interpolations.R2

    @time R_factor = compute_R_factor_current(trunk_sim.trunk, sol)
    plot_outer_trunk!(
        trunk,
        sol,
        R_factor;
        n = n,
        flipped = flipped,
        opacity = opacity_trunk
    )
    for m in muscle_indices
        opacity_f = opacity_fibers
        transp = false

        if m[3] == 1
            plot_muscle!(
                sol,
                R1[m[2]](m[1]),
                R2[m[2]](m[1]),
                trunk.Θ1R[m[1], m[2]],
                trunk.Θ2R[m[1], m[2]],
                trunk.Z1[m[1]],
                trunk.Z2[m[1]],
                colors[m[1], m[2], m[3]],
                R_factor;
                n = n,
                R1_surface = true,
                Z1_cap = (m[1] == 1),
                flipped = flipped,
                opacity_fibers = opacity_f,
                transparency = transp,
                override_transparent = false
            )
        else
            plot_muscle!(
                sol,
                R1[m[2]](m[1]),
                R2[m[2]](m[1]),
                trunk.Θ1L[m[1], m[2]],
                trunk.Θ2L[m[1], m[2]],
                trunk.Z1[m[1]],
                trunk.Z2[m[1]],
                colors[m[1], m[2], m[3]],
                R_factor;
                n = n,
                R1_surface = true,
                Z1_cap = (m[1] == 1),
                flipped = flipped,
                opacity_fibers = opacity_f,
                transparency = transp,
                override_transparent = false
            )
        end
    end
end

function plot_muscle_exploded!(
        R1,
        R2,
        Θ1,
        Θ2,
        Z1,
        Z2,
        offset_factor,
        color;
        n = 40,
        flipped = false,
        opacity = 1.0
)
    Θ = LinRange(Θ1, Θ2, n)
    Z = LinRange(Z1, Z2, n)
    ex = [1.0, 0.0, 0.0]
    offset = offset_factor * R2(Z1) * AngleAxis((Θ1 + Θ2) / 2.0, 0.0, 0.0, 1.0) * ex

    # R2 surface
    out = [[0.0, 0.0, u] + R2(u) * AngleAxis(v, 0.0, 0.0, 1.0) * ex + offset
           for u in Z,
    v in Θ]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = true,
        specular = Vec3f(0.6),
        transparency = true
    )

    # R1 surface
    out = [[0.0, 0.0, u] + R1(u) * AngleAxis(v, 0.0, 0.0, 1.0) * ex + offset
           for u in Z,
    v in Θ]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = false,
        specular = Vec3f(0.6),
        transparency = true
    )

    # Z1 cap
    R = LinRange(R1(Z1), R2(Z1), n)
    out = [[0.0, 0.0, Z1] + v * AngleAxis(u, 0.0, 0.0, 1.0) * ex + offset
           for u in Θ, v in R]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = false,
        specular = Vec3f(0.6),
        transparency = true
    )

    # Z2 cap
    R = LinRange(R1(Z2), R2(Z2), n)
    out = [[0.0, 0.0, Z2] + v * AngleAxis(u, 0.0, 0.0, 1.0) * ex + offset
           for u in Θ, v in R]
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = true,
        specular = Vec3f(0.6),
        transparency = true
    )

    # Θ1 side
    out = Matrix{SVector{3, Float64}}(undef, n, n)
    for i in eachindex(Z)
        u = Z[i]
        out[i, :] = [[0.0, 0.0, u] + v * AngleAxis(Θ1, 0.0, 0.0, 1.0) * ex + offset
                     for
                     v in LinRange(R1(u), R2(u), n)]
    end
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = true,
        specular = Vec3f(0.6),
        transparency = true
    )

    # Θ2 side
    out = Matrix{SVector{3, Float64}}(undef, n, n)
    for i in eachindex(Z)
        u = Z[i]
        out[i, :] = [[0.0, 0.0, u] + v * AngleAxis(Θ2, 0.0, 0.0, 1.0) * ex + offset
                     for
                     v in LinRange(R1(u), R2(u), n)]
    end
    x, y, z = extract_xyz(out, flipped)
    GLMakie.surface!(
        x,
        y,
        z,
        color = fill((color, opacity), n, n),
        ssao = true,
        invert_normals = false,
        specular = Vec3f(0.6),
        transparency = true
    )
end

"""
    $(TYPEDSIGNATURES)

Plots an exploded view of the muscles of the trunk `trunk_sim` for
the activated trunk quantities `a` of type `ActivatedTrunkQuantities`,
and the maximum activation magnitude `max_abs_γ` across all muscles. 
"""
function ActiveFilaments.plot_trunk_exploded!(trunk_sim::TrunkFast{T, N},
        a::ActivatedTrunkQuantities{T, N}, max_abs_γ::Float64; n = 40,
        colors = [:orange, :magenta, :cyan, :red, :blue],
        inactive_color = RGB(0.6, 0.6, 0.6), offset_factors = [1.5, 1.5, 1.2, 0.8, 0.5],
        flipped = false) where {T, N}
    trunk = trunk_sim.trunk
    R1 = trunk_sim.interpolations.R1
    R2 = trunk_sim.interpolations.R2
    γ = a.γ

    for i in 1:T
        for j in 1:5
            γR = γ[1][i, j]
            color_R = blend_color(inactive_color, colors[j], abs(γR) / max_abs_γ)
            opacity_R = 0.3 + abs(γR) / max_abs_γ * 0.7
            plot_muscle_exploded!(
                R1[j](i),
                R2[j](i),
                trunk.Θ1R[i, j],
                trunk.Θ2R[i, j],
                trunk.Z1[i],
                trunk.Z2[i],
                offset_factors[j],
                color_R;
                n = n,
                flipped = flipped,
                opacity = opacity_R
            )

            γL = γ[2][i, j]
            color_L = blend_color(inactive_color, colors[j], abs(γL) / max_abs_γ)
            opacity_L = 0.3 + abs(γL) / max_abs_γ * 0.7
            plot_muscle_exploded!(
                R1[j](i),
                R2[j](i),
                trunk.Θ1L[i, j],
                trunk.Θ2L[i, j],
                trunk.Z1[i],
                trunk.Z2[i],
                offset_factors[j],
                color_L;
                n = n,
                flipped = flipped,
                opacity = opacity_L
            )
        end
    end
end

function blend_color(c1, c2::Symbol, t)
    c2 = Colors.parse(Colorant, c2)
    return RGB(c1.r + t * (c2.r - c1.r), c1.g + t * (c2.g - c1.g), c1.b + t * (c2.b - c1.b))
end

end