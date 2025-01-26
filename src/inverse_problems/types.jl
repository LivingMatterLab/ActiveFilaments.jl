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

struct ConfigurationControlObjective{P<:DataType} <: AbstractControlObjective
    args::Vector{<:SArray}
    properties::Union{<:Any}
    propertyTypes::Vector{P}
    weights::Vector{<:SArray}
    ConfigurationControlObjective(
        args::Vector{<:SArray},
        properties::Vector{<:Any},
        weights::Vector{<:SArray};
        propertyTypes::Vector{P} = Vector{DataType}([r]),
    ) where {P<:DataType} =
        new{P}(args, properties, propertyTypes, weights)
end