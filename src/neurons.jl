include("calciumDynamics.jl")
include("ionChannels.jl")

@component function NaKaNeuron(;name, input_current=0.0, kwargs...)

    pars = @parameters begin
        C = 0.01
    end

    vars = @variables begin
        v(t) = -70
    end

    systems = @named begin
        na = HHSodiumChannel(v_in=v, parent=nothing)
        ka = ATypePotassiumChannel(v_in=v, parent=nothing)
    end

    eqs = [
        C * D(v) ~ na.base.i + ka.base.i
    ]

    return ODESystem(eqs, t, vars, pars; systems=systems, name=name)
end


@component function HHNeuron(;name, input_current=0.0, kwargs...)

    pars = @parameters begin
        C = 0.01
    end

    vars = @variables begin
        v(t) = -70
        #ca_influx() = 0.0
        #i(t), Will be useful and necessary when modeling networks, synapses
    end

    @named ca_dynamics = CalciumDynamics()

    systems = @named begin
        na = HHSodiumChannel(v_in=v, conductance=100.0, reversal_potential=50.0, parent=ca_dynamics)
        cas = SlowCalciumChannel(v_in=v, parent=ca_dynamics.E_Ca)
        cat = TransientCalciumChannel(v_in=v, conductance=1.3, parent=ca_dynamics.E_Ca)
        ka = ATypePotassiumChannel(v_in=v, conductance=5.0, reversal_potential=-80.0, parent=ca_dynamics)
        kca = CalciumActivatedPotassiumChannel(v_in=v, parent=ca_dynamics.Ca)
        kdr = DelayedRectifierPotassiumChannel(v_in=v, conductance=20.0, reversal_potential=-80.0, parent=ca_dynamics)
        h = HCurrentChannel(v_in=v, conductance=0.5, reversal_potential=-20.0, parent=ca_dynamics)
        #leak = LeakChannel(v_in=v, conductance=0.01, reversal_potential=-50.0, parent=ca_dynamics)
    end
    
    push!(systems, ca_dynamics)
    eqs = [
        C * D(v) ~ summ_channels(systems) + input_current,
        D(ca_dynamics.Ca) ~ (1/ca_dynamics.tau_Ca) * (ca_dynamics.Ca_base + 0.94*(cas.base.i + cat.base.i) - ca_dynamics.Ca)
    ]
    return ODESystem(eqs, t, vars, pars; systems=systems, name=name)
end

function summ_channels(systems)
    summ = sum(entity.base.i for entity in systems if hasproperty(entity, :base))
    return summ
end
function flux(channel)
    if hasproperty(channel, :I_Ca)
        return channel.I_Ca
    end
    return 0
end

#Add IonChannelFactory: registers instances of dynamicionchannel to an array at the neuron level, neuron sums it,
#Feeds that straight into it's CalciumDynamics
#Also removes the need manually write summing equation per neuron.