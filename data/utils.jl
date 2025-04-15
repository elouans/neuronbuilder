module utils

export load_params, get_unique_name


const NAME_DICT = Dict{DataType, Int}()
function get_unique_name(obj)
    T = typeof(obj)
    if !haskey(NAME_DICT, T)
        NAME_DICT[T] = 0
    end
    NAME_DICT[T] += 1
    return Symbol(lowercase(string(nameof(T))), "_", NAME_DICT[T])
end

#To be replaced with: load_specific_param_subset, etc

function load_params(filename::String="params.toml")
    param_path = joinpath(@__DIR__, "..", filename)
    return TOML.parsefile(param_path)
end
end