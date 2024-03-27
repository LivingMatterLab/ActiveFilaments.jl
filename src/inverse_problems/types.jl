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

# struct ConfigurationControlObjective{T, P <: ConfigurationProperty} <: AbstractControlObjective
#     args::SVector{T, Float64}
#     properties::Union{SVector{T, Float64}, SMatrix{T, 3, Float64}}
#     propertyType::Type{P}
#     ConfigurationControlObjective(args::SVector{T, Float64}, properties::SVector{T, Float64}, propertyType::Type{P}) where {T, P <: ScalarConfigurationProperty} = new{T, P}(args, properties, propertyType);
#     ConfigurationControlObjective(args::SVector{T, Float64}, properties::SMatrix{T, 3, Float64}; propertyType::Type{P} = r) where {T, P <: Vector3ConfigurationProperty} = new{T, P}(args, properties, propertyType);
# end

# struct ConfigurationControlObjective{T, S<:ConfigurationProperty, P <: Vector{S}} <: AbstractControlObjective
#     args::Vector{SVector{T, Float64}}
#     properties::Union{Vector{SVector{T, Float64}}, Vector{SMatrix{T, 3, Float64}}}
#     propertyTypes::Vector{S}
#     ConfigurationControlObjective(
#         args::Vector{SVector{T, Float64}}, 
#         properties::Vector{SVector{T, Float64}}, 
#         propertyTypes::Vector{S}
#         ) where {T, S <: ScalarConfigurationProperty} = 
#             new{T, S, Vector{S}}(args, properties, propertyTypes);
#     ConfigurationControlObjective(
#         args::Vector{SVector{T, Float64}}, 
#         properties::Vector{SMatrix{T, 3, Float64, U}}; 
#         propertyTypes::Vector{S} = Vector{Vector3ConfigurationProperty}([r])
#         ) where {T, S <: Vector3ConfigurationProperty, U} = 
#             new{T, S, Vector{S}}(args, properties, propertyTypes);
# end

struct ConfigurationControlObjective{P<:DataType} <: AbstractControlObjective
    args::Vector{<:SArray}
    properties::Union{<:Any}
    propertyTypes::Vector{P}
    weights::Vector{<:SArray}
    ConfigurationControlObjective(
        args::Vector{<:SArray}, 
        properties::Vector{<:Any},
        weights::Vector{<:SArray}; 
        propertyTypes::Vector{P} = Vector{DataType}([r])
        ) where {P <: DataType} = 
            new{P}(args, properties, propertyTypes, weights);
end