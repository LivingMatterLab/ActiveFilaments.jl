module PlottingExt
using ActiveFilaments
using ColorSchemes
using CairoMakie
using GLMakie

export 
    plotReachabilityCloudRGB

function plotReachabilityCloudRGB(sols, activationsGamma::Matrix{Float64}, gammaBounds, axesLimits; 
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

    x = getindex.(points, 1) * (flipped ? -1 : 1);
    y = getindex.(points, 2);
    z = getindex.(points, 3) * (flipped ? -1 : 1);

    activations = [activationsGamma[i, :] for i in axes(activationsGamma, 1)];

    gammaBoundsPairs = Vector{Tuple{Float32, Float32}}();
    for j in eachindex(gammaBounds)
        for k in eachindex(gammaBounds[j][1])
            bounds = Tuple(sort(collect(abs.([gammaBounds[j][1][k], gammaBounds[j][2][k]]))));
            push!(gammaBoundsPairs, bounds);
        end
    end
    # resX = 3840; resY = 2160;
    GLMakie.activate!()
    fig = Figure(resolution = resolution)
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = azimuth, elevation = elevation, perspectiveness = perspectiveness; kwargs...)
    GLMakie.scatter!(x, y, z, markersize = markersize, 
                    color = [
                        RGBA(activations[1][i] / gammaBoundsPairs[1][2], 
                        activations[2][i] / gammaBoundsPairs[2][2], 
                        activations[3][i] / gammaBoundsPairs[3][2], 
                        opacity) 
                        for i in 1:length(activations[1])
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

end