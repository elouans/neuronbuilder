module full_hh_single


using MTKNeuralComponents: HHNeuron
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using ModelingToolkitStandardLibrary.Blocks: Constant, RealInput
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots

export full_HH

function full_HH(input_current=1.0, plotting=true)
    println(input_current)
    @parameters t
    Neurone = HHNeuron(name=:neur, input_current=input_current)
    Neurone = structural_simplify(Neurone)

    prob = ODEProblem(Neurone, [], (0.0, 10))
    sol = solve(prob, Tsit5(), maxiters=1e9)

    println("First few timepoints: ", sol.t[1:5])
    println("First few voltage values: ", sol[Neurone.v][1:5])

    println("Something polite")
    if plotting
        println("Jean Plot Van Debugamme")
        p = plot(layout=(7,1), size=(1200,1200))

        # Plot membrane potential on top panel
        plot!(p[1], sol, vars=[Neurone.v], 
            label="Membrane Potential", 
            xlabel="Time (ms)", 
            ylabel="Voltage (mV)",
            linewidth=2,
            title="Membrane Potential")
        
        # Plot all currents on bottom panel
        plot!(p[2], sol, vars=[Neurone.ka.base.i,Neurone.h.base.i, Neurone.kdr.base.i], 
            label=["KCa" "H" "Kdr"], 
            xlabel="Time (ms)", 
            ylabel="Current",
            linewidth=2,
            title="Ionic Currents")

        plot!(p[3], sol, vars=[Neurone.cas.base.i, Neurone.cat.base.i], 
        label=["Cas_i" "cat_i"], 
        xlabel="Time (ms)", 
        ylabel="Current",
        linewidth=2,
        title="Calcium Channel I monitoring")

        plot!(p[4], sol, vars=[Neurone.ca_dynamics.E_Ca], 
        label=["ECA"], 
        xlabel="Time (ms)", 
        ylabel="Current",
        linewidth=2,
        title="Calcium Pot")

        plot!(p[5], sol, vars=[Neurone.ca_dynamics.Ca], 
        label=["Ca"], 
        xlabel="Time (ms)", 
        ylabel="Current",
        linewidth=2,
        title="Calcium Concentration")
        
        plot!(p[6], sol, vars=[Neurone.ka.base.i, Neurone.ka.base.m, Neurone.ka.base.h], 
        label=["kai" "kam" "kah"], 
        xlabel="Time (ms)", 
        ylabel="Current",
        linewidth=2,
        title="K Monitoring")

        plot!(p[7], sol, vars=[Neurone.na.base.i], 
        label=["Na"], 
        xlabel="Time (ms)", 
        ylabel="Current",
        linewidth=2,
        title="Na Monitoring")
        # Save to file
        savefig(p, "neuron_plot.png")
        
        # Optionally display in terminal if using Plots with GR backend
        display(p)
    end
    println("Hopin for a dream")
end
end