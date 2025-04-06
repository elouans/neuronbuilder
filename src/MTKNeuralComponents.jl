module MTKNeuralComponents

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq

include("../data/utils.jl")
include("calciumDynamics.jl")
include("ionChannels.jl")
include("neurons.jl")
include("../scripts/full_hh_single.jl")


export HHNeuron, dostuff

end 
