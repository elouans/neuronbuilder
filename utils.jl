module NeuronUtils

export load_params, get_unique_name, calc_ca_reversal

using TOML

"""
    get_unique_name(obj)

Generate a unique name for an object based on its type.
Maintains a counter for each object type to ensure uniqueness.
"""
const NAME_DICT = Dict{DataType, Int}()
function get_unique_name(obj)
    T = typeof(obj)
    if !haskey(NAME_DICT, T)
        NAME_DICT[T] = 0
    end
    NAME_DICT[T] += 1
    return Symbol(lowercase(string(nameof(T))), "_", NAME_DICT[T])
end

"""
    load_params(filename::String)

Load model parameters from a TOML configuration file.
"""
function load_params(filename::String)
    return TOML.parsefile(filename)
end
end