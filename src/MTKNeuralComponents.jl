module MTKNeuralComponents

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq

include("utils.jl")

include("neurons.jl")

export HHNeuron, NaKaNeuron

end 
