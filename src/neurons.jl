module neurons

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using ..calciumDynamics
using ..ionChannels

export HHNeuron


@component function HHNeuron(;name, input_current=0.0, kwargs...)

    pars = @parameters begin
        C = 0.01
        gap_C = 0.02
        #spike_threshold=-40.0 #Useless now, artifact from modeling spikes as discrete events
    end

    vars = @variables begin
        v(t) = -70
        #i(t), Will be useful and necessary when modeling networks, synapses
    end

    @named ca_dynamics = calciumDynamics.CalciumDynamics(v_in=v)

    systems = @named begin
        na = ionChannels.HHSodiumChannel(v_in=v, conductance=100.0, reversal_potential=50.0)
        cas = ionChannels.SlowCalciumChannel(v_in=v, conductance=3.0, ca_dynamics=ca_dynamics)
        cat = ionChannels.TransientCalciumChannel(v_in=v, condutance=1.3, ca_dynamics=ca_dynamics)
        ka = ionChannels.ATypePotassiumChannel(v_in=v, conductance=5.0, reversal_potential=-80.0)
        kca = ionChannels.CalciumActivatedPotassiumChannel(v_in=v, ca_dynamics=ca_dynamics)
        kdr = ionChannels.DelayedRectifierPotassiumChannel(v_in=v, conductance=20.0, reversal_potential=-80.0)
        h = ionChannels.HCurrentChannel(v_in=v, conductance=0.5, reversal_potential=-20.0)
        leak = ionChannels.LeakChannel(v_in=v, condutance=0.01, reversal_potential=-50.0)
    end
    
    #push!(systems, ca_dynamics)
    #Simplify, sum all channels' intensities
    eqs = [
        C * D(v) ~ na.base.i + cas.base.pin.i + cat.base.pin.i + ka.base.i + kdr.base.i + h.base.i + leak.base.i + input_current
        #C * D(v) ~ na.base.i+ ka.base.i + kdr.base.i + h.base.i + leak.base.i + input_current

    ]
    return ODESystem(eqs, t, vars, pars; systems=systems, name=name)
end
end