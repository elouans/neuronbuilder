module neurons

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using ..calciumDynamics
using ..utils
using ..ionChannels

export HHNeuron


@component function HHNeuron(;name, input_current=0.0, kwargs...)
    params = utils.load_params()

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
        na = ionChannels.HHSodiumChannel(v_in=v, conductance=params.gNabar, reversal_potential=params.eNa)
        cas = ionChannels.SlowCalciumChannel(v_in=v, conductance=params.gCaSbar, ca_dynamics=ca_dynamics)
        cat = ionChannels.TransientCalciumChannel(v_in=v, condutance=params.gCaTbar, ca_dynamics=ca_dynamics)
        ka = ionChannels.ATypePotassiumChannel(v_in=v, conductance=params.gKabar, reversal_potential=params.eK)
        #kca = ionChannels.CalciumActivatedPotassiumChannel(v_in=v, ca_dynamics=ca_dynamics)
        kdr = ionChannels.DelayedRectifierPotassiumChannel(v_in=v, conductance=params.gKdrbar, reversal_potential=eK)
        h = ionChannels.HCurrentChannel(v_in=v, conductance=gHbar, reversal_potential=eH)
        leak = ionChannels.LeakChannel(v_in=v, condutance=gleak, reversal_potential=eleak)
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