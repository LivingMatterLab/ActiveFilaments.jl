#########################################################################
### Abstract supertypes for the Active Filament theory implementation ###
#########################################################################
#region ===========================
"""
Abstract supertype for various filament design properties
"""
abstract type AbstractDesignProperty end

"""
Abstract supertype for the mechanical architecture of the filament
"""
abstract type AbstractMechanicalProperties <: AbstractDesignProperty end

"""
Abstract supertype for filament geometry
"""
abstract type AbstractGeometry <: AbstractDesignProperty end

"""
Abstract supertype for fiber architecture of the filament
"""
abstract type AbstractFiberArchitecture <: AbstractDesignProperty end

"""
Abstract supertype for rings
"""
abstract type AbstractRing end

"""
Abstract supertype for the activation pattern
"""
abstract type AbstractActivation end

"""
Abstract supertype for the piecewise activation pattern
"""
abstract type AbstractActivationPiecewise <: AbstractActivation end

"""
Abstract supertype for the precomputed quantities structure
"""
abstract type AbstractPrecomputedQuantities end
#endregion ===========================