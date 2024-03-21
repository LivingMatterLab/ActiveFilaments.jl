"""
Main module based on the Active Filament theory.

Packages used in the module:
$(IMPORTS)

# Exports:
$(EXPORTS)
"""
module ActiveFilaments

######### Imported packages
#region ===========================
using Preferences, UUIDs
using Revise;
using Parameters;
using StaticArrays;
#Preferences.set_preferences!(UUID("764a87c0-6b3e-53db-9096-fe964310641d"),
#    "PrecompileShooting" => false)
using DifferentialEquations;
# using BoundaryValueDiffEq;
using Plots;
using Rotations;
using LinearAlgebra;
using BenchmarkTools;
using CSV;
using Tables;
using DataFrames;
using JLD2;
using Symbolics;
using SymbolicNumericIntegration;
# using DiffEqGPU;
using Distributions;
using DocStringExtensions;
using QuadGK;
using Interpolations;
using Optimization;
using OptimizationOptimJL;
using Distances;
using ForwardDiff;
# using Zygote;
using ReverseDiff;
# using ModelingToolkit;
using OptimizationNLopt;

import Base.copy
#endregion ========================

######### Module exports
export 
    PiecewiseStructure,
    MechanicalProperties,
    Geometry,
    FiberArchitecture,
    Ring,
    AuxiliaryProperties,
    AFilament,
    InnerTube,
    FilamentStiffness,
    ActivationPiecewiseGamma,
    ActivationFourier,
    FiberID,

    piecewiseGammaToFourier,
    computePropertyPrefactors,
    computeUQuantities,
    computeUHat,
    generatePrecomputedQuantitiesSA,
    convertUQuantToStatic,

    selfWeightSolve,
    selfWeightSolveSym,
    solveIntrinsic,
    generateIntrinsicReachVol,
    generateSelfWeightReachVol,
    generateSelfWeightReachVolSym,

    loopRadius,
    bendingAngle,
    tiltAngle,
    tiltAngleCapped,
    fiberStrain,

    ConfigurationControlObjective,
    optimizeActivation,
    r,
    d1,
    d2,
    d3,
    build_distance_function,

    Trunk,
    TrunkFast,
    ActivatedTrunkQuantities,
    constant_activation,
    computeH,
    self_weight_solve,
    build_trunk_bvp,
    self_weight_solve_single,
    optimize_activation,
    rotate_bc,

    plotReachabilityCloudRGB,
    plotReachabilityCloud,
    plotReachabilityCloudRGBSlice,
    plotReachabilityCloudSlice,
    plotFilamentCollapsedRings!,
    plotConfigurationsSelfWeight!,
    plotConfigurationTubesSelfWeight!,

    plot_trunk!,
    plot_trunk_exploded!

# Main scripts
include("types/types.jl");
include("model/model.jl");
include("trunk/trunk.jl")
include("computation/computation.jl");
include("utils.jl");
include("analysis/analysis.jl");
include("inverse_problems/inverse.jl")
include("legacy/legacy.jl") # Legacy (comment out if unnecessary): needed for some legacy data imports

# Plotting functionality
include("plotting/plotting.jl")
    

end