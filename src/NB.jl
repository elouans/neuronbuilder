module NB

include("../data/utils.jl")
include("calciumDynamics.jl")
include("ionChannels.jl")
include("neurons.jl")
include("../scripts/dbg_1n.jl")

# Re-export key components
using .utils
using .calciumDynamics
using .ionChannels
using .neurons

export HHNeuron, neurons

end # module NB
