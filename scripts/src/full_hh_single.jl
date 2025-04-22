module full_hh_single


using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using MTKNeuralComponents: HHNeuron
using Plots

export full_HH

test() = print("Hello World!")
function full_HH(input_current=10.0, plotting=true)
    println(input_current)
    @parameters t
    Neurone = HHNeuron(name=:neur, input_current=input_current)
    println(size(ModelingToolkit.unknowns(Neurone)))
    println(Neurone.discrete_subsystems)
    Neurone = structural_simplify(Neurone)

    prob = ODEProblem(Neurone, [], (0.0, 40))
    sol = solve(prob, Tsit5())
    println("Something polite")
    if plotting
        println("Jean Plot Van Debugamme")
        p = plot(layout=(3,1), size=(800,600))

        # Plot membrane potential on top panel
        plot!(p[1], sol, vars=[Neurone.v], 
            label="Membrane Potential", 
            xlabel="Time (ms)", 
            ylabel="Voltage (mV)",
            linewidth=2,
            title="Membrane Potential")

        # Plot all currents on bottom panel
        plot!(p[2], sol, vars=[Neurone.na.base.i, Neurone.cas.base.i, Neurone.cat.base.i, 
                            Neurone.ka.base.i, Neurone.kca.base.i, Neurone.leak.base.i], 
            label=["Na+" "CaS" "CaT" "KA" "KCa" "Leak"], 
            xlabel="Time (ms)", 
            ylabel="Current",
            linewidth=2,
            title="Ionic Currents")

        plot!(p[3], sol, vars=[Neurone.ca_dynamics.I_Ca, Neurone.ca_dynamics.E_Ca, Neurone.ca_dynamics.Ca], 
        label=["ICA" "ECA" "Ca"], 
        xlabel="Time (ms)", 
        ylabel="Current",
        linewidth=2,
        title="Calcium Dynamics")
        
        # Save to file
        savefig(p, "neuron_plot.png")
        
        # Optionally display in terminal if using Plots with GR backend
        display(p)
    end
    println("Hopin for a dream")
end
end