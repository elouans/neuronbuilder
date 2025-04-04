module NB

include("../data/utils.jl")
include("calciumDynamics.jl")
include("ionChannels.jl")
include("neurons.jl")
include("../scripts/full_hh_single.jl")

# Re-export key components
using .utils
using .calciumDynamics
using .ionChannels
using .neurons

export HHNeuron, dostuff

end # module NB
