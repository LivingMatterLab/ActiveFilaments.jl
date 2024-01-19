# ActiveFilaments.jl
ActiveFilaments.jl is a computationally efficient Julia implementation of the mechanical theory of active slender structures found in:

[Kaczmarski, B., Moulton, D. E., Kuhl, E. & Goriely, A. Active filaments I: Curvature and torsion generation. _Journal of the Mechanics and Physics of Solids_ **164**, 104918 (2022).](https://doi.org/10.1016/j.jmps.2022.104918)

Using this package, you can design the geometry and fiber architecture in a slender structure and simulate its deformation due to fibrillar activation. The package provides a library of analysis tools for exploratory soft-robotics research including fast workspace computation for active soft slender structures.

# Installation
Run the following command in Julia 1.10.0+:
```julia
] add https://github.com/LivingMatterLab/ActiveFilaments.git
```
Wait for all dependencies to precompile. Slowdowns are possible in parallelized BVP solutions if Julia version 1.9.0 or older is used.

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
```
