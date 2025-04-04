module full_hh_single


using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using NB
using Plots

export dostuff

test() = print("Hello World!")
function dostuff()
    Neurone = NB.neurons.HHNeuron(name=:neur)
    println(size(ModelingToolkit.unknowns(Neurone)))
    println(Neurone.discrete_subsystems)
    Neurone = structural_simplify(Neurone)

    prob = ODEProblem(Neurone, [], (0.0, 20))
    sol = solve(prob, Tsit5())
    print(ln("SHIT"))
end
end