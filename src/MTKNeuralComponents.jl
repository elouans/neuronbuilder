module MTKNeuralComponents

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant, RealInput
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkit: get_eqs
using OrdinaryDiffEq

include("utils.jl")

include("neurons.jl")


export HHNeuron, full_HH, minimal_HH, NaKaNeuron

end 
