# ######### Fiber activation structures
#region ===========================
# Might change Num to Float32 below to check for compatibility with CUDA
"""
    $(TYPEDEF)

Piecewise fiber activation structure (for GPU computation).

$(TYPEDFIELDS)
"""
mutable struct ActivationPiecewiseGammaGPU <: AbstractActivationPiecewise
    "Number of angular sectors"
    N::Int
    "Array of activation magnitudes for each sector"
    γ::Vector{Num}
    "Angular extent of each sector"
    σ::Num
    "Angular offset of each sector"
    θ0::Num
end

"""
    $(TYPEDEF)

Piecewise fiber activation structure (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
mutable struct ActivationPiecewiseGamma <: AbstractActivationPiecewise
    "Number of angular sectors"
    N::Int
    "Array of activation magnitudes for each sector"
    γ
    "Angular extent of each sector"
    σ
    "Angular offset of each sector"
    θ0
end

"""
    $(TYPEDEF)

Fourier coefficients of fiber activation (for GPU computation).

$(TYPEDFIELDS)
"""
struct ActivationFourierGPU <: AbstractActivationPiecewise
    "Coefficient 1"
    a0::Num
    "Coefficient 2"
    a1::Num
    "Coefficient 3"
    b1::Num
end

"""
    $(TYPEDEF)

Fourier coefficients of fiber activation (symbolic or numeric computation).

$(TYPEDFIELDS)
"""
struct ActivationFourier{T} <: AbstractActivationPiecewise
    # Assumes constant a0 (this is the case if the number of fibers, 
    # σ etc. doesn't change with Z); a0 is still constant if α2 or θ0 changes
    "Coefficient 1"
    a0::Float64 
    "Coefficient 2"
    a1::T
    "Coefficient 3"
    b1::T
end
#endregion ===========================
