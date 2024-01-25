# ActiveFilaments.jl
ActiveFilaments.jl is a computationally efficient Julia implementation of the mechanical theory of active slender structures found in:

[Kaczmarski, B., Moulton, D. E., Kuhl, E. & Goriely, A. Active filaments I: Curvature and torsion generation. _Journal of the Mechanics and Physics of Solids_ **164**, 104918 (2022).](https://doi.org/10.1016/j.jmps.2022.104918)

Using this package, you can design the geometry and fiber architecture in a slender structure and simulate its deformation due to fibrillar activation. The package provides a library of analysis tools for exploratory biomechanics and soft-robotics research including fast workspace computation for active soft slender structures.

# Installation
Run the following command in Julia 1.10.0+:
```
] add https://github.com/LivingMatterLab/ActiveFilaments.git
```
Wait for all dependencies to precompile. Slowdowns are possible during precompilation in parallelized BVP solutions if Julia version 1.9.0 or older is used.

# Example usage
In the simplest case of uniform material properties and a single fiber ring, we can define an active filament by running
```julia
using ActiveFilaments

# Mechanical properties of the ring
E = 1.0e6 # Young's modulus
ν = 0.5   # Incompressible
ρ = 1000.0 # Volumetric density
mech = MechanicalProperties(E, ν)

# Geometry of the ring
L = 1.0 # Length of the filament
R2 = L / 20.0 # Outer radius of the ring
R1 = R2 * 0.8 # Inner radius of the ring
geom = Geometry(R1, R2)

# Fiber architecture of the ring
α2 = pi / 8.0 # Outer helical angle of the fiber architecture
arch = FiberArchitecture(α2)

# Define the filament structure
rings = [Ring(mech, geom, arch)]
filament = AFilament(rings; L = L, ρvol = ρ)
```
We can then prescribe a piecewise constant activation pattern with three independent 30° activation sectors uniformly spaced around the ring starting from 0° at the +_x_-axis. Let us choose contractile fiber activations (-2.0, -1.3, 0.0) for the three respective sectors: 
```julia
activation = [ActivationPiecewiseGamma(3, [-2.0, -1.3, 0.0], 30.0 / 180.0 * pi, 0.0)]  
```
The deformation due to the above activation without external loading can then be computed as
```julia
sol = solveIntrinsic(filament, activation)
```
The resulting `sol` structure contains the following interpolating functions:
$$\texttt{sol}(Z) = \big[r_X(Z), r_Y(Z), r_Z(Z), d_{1X}(Z), d_{1Y}(Z), d_{1Z}(Z), d_{2X}(Z), d_{2Y}(Z), d_{2Z}(Z), d_{3X}(Z), d_{3Y}(Z), d_{3Z}(Z)\big]$$
For example, evaluating `sol` at $L / 10$ gives
```
julia> sol(L / 10.0)
12-element StaticArraysCore.SVector{12, Float64} with indices SOneTo(12):
  0.008311388337878088
  0.009719221098605893
  0.09072971511373011
  0.9421455745518371
 -0.32511595347415
 -0.0816169905655878
  0.294271010395179
  0.9188022897556413
 -0.2630720904634142
  0.16051881133807738
  0.2238346915390569
  0.9613177113058022
```
