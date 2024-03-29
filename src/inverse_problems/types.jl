"Abstract supertype for the control objective"
abstract type AbstractControlObjective end

"Abstract supertype for configuration properties"
abstract type ConfigurationProperty end

"Abstract supertype for 3D vector configuration properties"
abstract type Vector3ConfigurationProperty <: ConfigurationProperty end

"Abstract supertype for scalar configuration properties"
abstract type ScalarConfigurationProperty <: ConfigurationProperty end

struct r <: Vector3ConfigurationProperty end
struct d1 <: Vector3ConfigurationProperty end
struct d2 <: Vector3ConfigurationProperty end
struct d3 <: Vector3ConfigurationProperty end
struct ζ <: ScalarConfigurationProperty end
struct u1 <: ScalarConfigurationProperty end
struct u2 <: ScalarConfigurationProperty end
struct u3 <: ScalarConfigurationProperty end

struct ConfigurationControlObjective{T, P <: ConfigurationProperty} <: AbstractControlObjective
    args::SVector{T, Float64}
    properties::Union{SVector{T, Float64}, SMatrix{T, 3, Float64}}
    propertyType::Type{P}
    ConfigurationControlObjective(args::SVector{T, Float64}, properties::SVector{T, Float64}, propertyType::Type{P}) where {T, P <: ScalarConfigurationProperty} = new{T, P}(args, properties, propertyType);
    ConfigurationControlObjective(args::SVector{T, Float64}, properties::SMatrix{T, 3, Float64}; propertyType::Type{P} = r) where {T, P <: Vector3ConfigurationProperty} = new{T, P}(args, properties, propertyType);
end

struct ConfigurationControlObjective{N, T, P <: SVector{N, <:ConfigurationProperty}} <: AbstractControlObjective
    args::SVector{N, SVector{T, Float64}}
    properties::Union{SVector{N, SVector{T, Float64}}, SVector{N, SMatrix{T, 3, Float64}}}
    propertyTypes::SVector{N, Type{P}}
    ConfigurationControlObjective(args::SVector{N, SVector{T, Float64}}, properties::SVector{N, SVector{T, Float64}}, propertyTypes::SVector{N, Type{P}}) where {N, T, P <: ScalarConfigurationProperty} = new{N, T, P}(args, properties, propertyTypes);
    ConfigurationControlObjective(args::SVector{N, SVector{T, Float64}}, properties::SVector{N, SMatrix{T, 3, Float64}}; propertyTypes::SVector{N, Type{P}} = SVector{1, Vector3ConfigurationProperty}([r])) where {N, T, P <: Vector3ConfigurationProperty} = new{N, T, P}(args, properties, propertyTypes);
end