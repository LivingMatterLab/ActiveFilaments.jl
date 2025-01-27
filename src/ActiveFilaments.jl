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
using Revise
using Parameters
using StaticArrays
using DifferentialEquations
using Plots
using Rotations
using LinearAlgebra
using BenchmarkTools
using CSV
using Tables
using DataFrames
using JLD2
using Symbolics
using SymbolicNumericIntegration
using Distributions
using DocStringExtensions
using QuadGK
using Interpolations
using Optimization
using Distances
using ForwardDiff
using OptimizationNLopt
using OptimizationBBO
using OptimizationNOMAD
using ODEInterfaceDiffEq

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
       FiberID, piecewiseGammaToFourier,
       computePropertyPrefactors,
       computeUQuantities,
       computeUHat,
       generatePrecomputedQuantitiesSA,
       convertUQuantToStatic, selfWeightSolve,
       selfWeightSolveSym,
       solveIntrinsic,
       generateIntrinsicReachVol,
       generateSelfWeightReachVol,
       generateSelfWeightReachVolSym, loopRadius,
       bendingAngle,
       tiltAngle,
       tiltAngleCapped,
       fiberStrain, ConfigurationControlObjective,
       optimizeActivation,
       r,
       d1,
       d2,
       d3,
       build_distance_function, Trunk,
       TrunkFast,
       ActivatedTrunkQuantities,
       constant_activation,
       computeH,
       self_weight_solve,
       build_trunk_bvp,
       build_trunk_wrap_bvp,
       build_trunk_wrap_nbodyc_bvp,
       build_trunk_wrap_offset_nbodyc_bvp,
       self_weight_solve_single,
       optimize_activation,
       rotate_bc,
       build_trunk_ivp,
       ivp_solve_single,
       compute_R_factor_current, plotReachabilityCloudRGB,
       plotReachabilityCloud,
       plotReachabilityCloudRGBSlice,
       plotReachabilityCloudSlice,
       plotFilamentCollapsedRings!,
       plotConfigurationsSelfWeight!,
       plotConfigurationTubesSelfWeight!, plot_trunk!,
       plot_trunk_exploded!,
       plot_trunk_isolated!

# Main scripts
include("types/types.jl")
include("model/model.jl")
include("trunk/trunk.jl")
include("computation/computation.jl")
include("utils.jl")
include("analysis/analysis.jl")
include("inverse_problems/inverse.jl")

# Legacy (comment out if unnecessary): needed for some legacy data imports
include("legacy/legacy.jl")

# Plotting functionality
include("plotting/plotting.jl")

end