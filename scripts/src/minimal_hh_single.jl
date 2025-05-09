module minimal_hh_single

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using MTKNeuralComponents: NaKaNeuron
using Plots

export minimal_HH

test() = println("Test")

function minimal_HH()
    Neurone = NaKaNeuron(name=:neur)
    Neurone = structural_simplify(Neurone)

    prob = ODEProblem(Neurone, [], (0.0, 20))
    sol = solve(prob, Tsit5())
    println("Something polite")

    p = plot(sol, vars=[Neurone.v], 
             label="Membrane Potential", 
             xlabel="Time (ms)", 
             ylabel="Voltage (mV)")
    
    # Save to file
    savefig(p, "neuron_plot.png")
    
    # Optionally display in terminal if using Plots with GR backend
    display(p)

    println("Hopin for a dream")
end
end