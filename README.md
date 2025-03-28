# ActiveFilaments.jl
ActiveFilaments.jl is a computationally efficient Julia implementation of the mechanical theory of active slender structures found in:

[Kaczmarski, B., Moulton, D. E., Kuhl, E. & Goriely, A. Active filaments I: Curvature and torsion generation. _Journal of the Mechanics and Physics of Solids_ **164**, 104918 (2022).](https://doi.org/10.1016/j.jmps.2022.104918)

By using this package, you can design the geometry and fiber architecture in a slender structure and simulate its deformation due to fibrillar activation. The package provides a library of analysis tools for exploratory biomechanics and soft-robotics research including fast workspace computation for active soft slender structures.

# Installation
Run the following command in Julia 1.10.0+:
```
] add https://github.com/LivingMatterLab/ActiveFilaments.git
```
Wait for all dependencies to precompile. Slowdowns are possible during precompilation in parallelized BVP solutions if Julia v1.9.0 or older is used.

# Example usage
## Computing an activated configuration
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
R2 = L / 16.0 # Outer radius of the ring
R1 = R2 * 0.6 # Inner radius of the ring
geom = Geometry(R1, R2)

# Fiber architecture of the ring
α2 = pi / 10.0 # Outer helical angle of the fiber architecture
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
$$\texttt{sol}(Z) = \big[r_X(Z), r_Y(Z), r_Z(Z), d_{1X}(Z), d_{1Y}(Z), d_{1Z}(Z), d_{2X}(Z), d_{2Y}(Z), d_{2Z}(Z), d_{3X}(Z), d_{3Y}(Z), d_{3Z}(Z)\big]$$.

For example, evaluating `sol` at $L / 10$ gives
```
julia> sol(L / 10.0)
12-element StaticArraysCore.SVector{12, Float64} with indices SOneTo(12):
  0.010887432435510243
  0.010588028536001756
  0.08230390930976722
  0.9378971336343721
 -0.3068269347479843
 -0.16188328769080979
  0.24481737934065784
  0.9160122985605073
 -0.31778281838178346
  0.24579141053584871
  0.25841575222701285
  0.9342418752698606
```
## Computing and plotting a reachability cloud
Let us continue with the filament defined in the previous example. For the three-sector design of the ring, we can define the bounds on the activation
sampling region as:
```julia
γBounds = [([-2.4, -2.0, -2.7], [0.0, 0.0, 0.0])]
```
Note that the array contains only one `Tuple` since the filament has a single fiber ring. In the above example, the activation bounds are such that
$$\gamma\in[-2.4,0]$$ in the first sector, $$\gamma\in[-2,0]$$, in the second sector, and $$\gamma\in[-2.7,0]$$ in the third sector.

We then specify the number of configurations to compute during cloud generation:
```julia
nTrajectories = 2000000
```
2 million configurations without external loading typically take from several seconds to a minute to compute, depending on the number of accessible threads.

We then indicate a custom directory and file path for the cloud data:
```julia
dir = "ReachabilityCloudOutput"
path = string(dir, "/one_ring_three_fibers")
```

We can now simply call `generateIntrinsicReachVol` to generate a reachability cloud for activated configurations (without external loading)
bounded by `γBounds`.
```julia
(sol, activationsGamma) = generateIntrinsicReachVol(
    filament, activation, γBounds, nTrajectories, path; 
    save_gamma_structs = false)
```
An output similar to the one shown below should appear (time and memory allocation values will be different):
```julia
Precomputation started...
  1.825260 seconds (190.00 M allocations: 7.615 GiB, 34.44% gc time)
Precomputation complete.
Computing cloud...
  2.811433 seconds (58.00 M allocations: 8.285 GiB)
Cloud computed!
Saving...
Saved!
```
In the `generateIntrinsicReachVol` function call, `activation` serves as an activation structure indicator for the reachability cloud generator. Setting the flag `save_gamma_structs` to `false` 
saves the activation samples in a simple `.csv` file.

To use the output data effectively, we can also save the filament and activation data created before:
```julia
using JLD2
@save string(dir, "/filament.jld2") filament
@save string(dir, "/gammaBounds.jld2") γBounds
@save string(dir, "/activationStructure.jld2") activation
```

In a new script, we can import all the data as follows:
```julia
using JLD2
using CSV
@load "ReachabilityCloudOutput/one_ring_three_fibers.jld2" sol
@load "ReachabilityCloudOutput/filament.jld2" filament
@load "ReachabilityCloudOutput/gammaBounds.jld2" γBounds
activationsGamma = CSV.read("ReachabilityCloudOutput/one_ring_three_fibers_gamma.csv" , 
                            CSV.Tables.matrix; header = false)
```

To visualize the cloud, we import `GLMakie` via `using GLMakie` and use the PlottingExt functionality of this package to generate the following animation:

https://github.com/user-attachments/assets/c71d4287-fd44-4bd8-b8bc-6c7986b2ffaa

In the above visualization, the points are color-coded using the `:thermal` colormap according to the activation magnitude in the third fiber.

# Other work that uses this package
[Kaczmarski, B., Moulton, D. E., Goriely, Z., Goriely, A. & Kuhl, E. Ultra-fast physics-based modeling of the elephant trunk. _bioRxiv_ 2024.10.27.620533 (2024).](https://doi.org/10.1101/2024.10.27.620533)

[Kaczmarski, B., Moulton, D. E., Goriely, A. & Kuhl, E. Minimal activation with maximal reach: Reachability clouds of bio-inspired slender manipulators. _Extreme Mechanics Letters_ **71**, 102207 (2024).](https://doi.org/10.1016/j.eml.2024.102207)

[Kaczmarski, B. et al. Minimal Design of the Elephant Trunk as an Active Filament. _Physical Review Letters_ **132**, 248402 (2024).](https://doi.org/10.1103/PhysRevLett.132.248402)
