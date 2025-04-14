module MTKNeuralComponents

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq

include("utils.jl")

include("neurons.jl")
include("../scripts/full_hh_single.jl")
include("../scripts/minimal_hh_single.jl")

function something()
    println("BBBB")
    full_HH()
    println("BBefefBB")
end

function check_files()
    println("Current directory: ", pwd())
    println("Does ../scripts/full_hh_single.jl exist? ", isfile("../scripts/full_hh_single.jl"))
    println("Does ../scripts/minimal_hh_single.jl exist? ", isfile("../scripts/minimal_hh_single.jl"))
end
export check_files

export HHNeuron, full_HH, minimal_HH, something

end 
