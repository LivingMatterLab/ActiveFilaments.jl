######################################################
### Types for the inverse problem solution methods ###
######################################################
"""
Abstract supertype for the control objective
"""
abstract type AbstractControlObjective end

"""
Abstract supertype for configuration properties
"""
abstract type ConfigurationProperty end

"""
Abstract supertype for 3D vector configuration properties
"""
abstract type Vector3ConfigurationProperty <: ConfigurationProperty end

"""
Abstract supertype for scalar configuration properties
"""
abstract type ScalarConfigurationProperty <: ConfigurationProperty end

"""
Defines the centerline property type.
"""
struct r <: Vector3ConfigurationProperty end
"""
Defines the d1 director property type.
"""
struct d1 <: Vector3ConfigurationProperty end
"""
Defines the d2 director property type.
"""
struct d2 <: Vector3ConfigurationProperty end
"""
Defines the d3 director property type.
"""
struct d3 <: Vector3ConfigurationProperty end
"""
Defines the extension property type.
"""
struct ζ <: ScalarConfigurationProperty end
"""
Defines the u1 curvature component property type.
"""
struct u1 <: ScalarConfigurationProperty end
"""
Defines the u2 curvature component property type.
"""
struct u2 <: ScalarConfigurationProperty end
"""
Defines the u3 curvature component property type.
"""
struct u3 <: ScalarConfigurationProperty end

"""
    $(TYPEDEF)

Defines the control objective in terms of the properties of the desired configuration.

$(TYPEDFIELDS)

Uses a constructor of the form: 
```
ConfigurationControlObjective(
            args::Vector{<:SArray},
            properties::Vector{<:Any},
            weights::Vector{<:SArray};
            propertyTypes::Vector{P} = Vector{DataType}([r])
) where {P <: DataType}
```
"""
struct ConfigurationControlObjective{P <: DataType} <: AbstractControlObjective
    "`Z` arguments at which the desired properties are defined"
    args::Vector{<:SArray}
    "The values of the properties at `args`"
    properties::Union{<:Any}
    "The property types associated with each element in `properties`"
    propertyTypes::Vector{P}
    "The weights associated with each property at each element of `args`"
    weights::Vector{<:SArray}
    function ConfigurationControlObjective(
            args::Vector{<:SArray},
            properties::Vector{<:Any},
            weights::Vector{<:SArray};
            propertyTypes::Vector{P} = Vector{DataType}([r])
    ) where {P <: DataType}
        new{P}(args, properties, propertyTypes, weights)
    end
end