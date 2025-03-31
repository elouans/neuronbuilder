module dbg_1n

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using NB
using Plots

test() = print("Hello World!")

Neurone = NB.neurons.HHNeuron(name=:neur)
println(size(ModelingToolkit.unknowns(Neurone)))
println(Neurone.discrete_subsystems)
Neurone = structural_simplify(Neurone)

prob = ODEProblem(Neurone, [], (0.0, 20))
sol = solve(prob, Tsit5())
end