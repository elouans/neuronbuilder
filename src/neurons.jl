include("calciumDynamics.jl")
include("ionChannels.jl")

@component function HHNeuron(;name, input_current=0.0, kwargs...)

    pars = @parameters begin
        C = 0.01
        gap_C = 0.02
        #spike_threshold=-40.0 #Useless now, artifact from modeling spikes as discrete events
    end

    vars = @variables begin
        v(t) = -70
        ca_influx() = 0.0
        #i(t), Will be useful and necessary when modeling networks, synapses
    end

    @named ca_dynamics = CalciumDynamics(v_in=v)

    systems = @named begin
        na = HHSodiumChannel(v_in=v, conductance=100.0, reversal_potential=50.0)
        cas = SlowCalciumChannel(v_in=v, conductance=3.0, ca_dynamics=ca_dynamics)
        cat = TransientCalciumChannel(v_in=v, conductance=1.3, ca_dynamics=ca_dynamics)
        ka = ATypePotassiumChannel(v_in=v, conductance=5.0, reversal_potential=-80.0)
        kca = CalciumActivatedPotassiumChannel(v_in=v, ca_dynamics=ca_dynamics)
        kdr = DelayedRectifierPotassiumChannel(v_in=v, conductance=20.0, reversal_potential=-80.0)
        h = HCurrentChannel(v_in=v, conductance=0.5, reversal_potential=-20.0)
        leak = LeakChannel(v_in=v, conductance=0.01, reversal_potential=-50.0)
    end
    
    #push!(systems, ca_dynamics)
    #Simplify, sum all channels' intensities
    eqs = [
        C * D(v) ~ na.base.i + cas.base.i + cat.base.i + ka.base.i  + kca.base.i + kdr.base.i + h.base.i + leak.base.i + input_current
        #C * D(v) ~ na.base.i+ ka.base.i + kdr.base.i + h.base.i + leak.base.i + input_current
        #ca_dynamics.I_Ca ~ cas.I_Ca + cat.I_Ca + kca.pin.i
        #C * D(v) ~ cat.base.i + input_current


    ]
    return ODESystem(eqs, t, vars, pars; systems=systems, name=name)
end
#Add IonChannelFactory: registers instances of dynamicionchannel to an array at the neuron level, neuron sums it,
#Feeds that straight into it's CalciumDynamics
#Also removes the need manually write summing equation per neuron.