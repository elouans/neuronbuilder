include("calciumDynamics.jl")
include("ionChannels.jl")


@mtkmodel InjectedCapacitor begin
    @components begin
        I = RealInput()
    end
    @parameters begin
        C, [description = "Capacitance"]
    end
    @extend v, i = oneport = OnePort(; v)
    @variables begin
        i_fake(t), [output=true]
    end
    @equations begin
        D(v) ~ (i + I.u) / C
        i_fake ~ i
    end
end

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
        C * D(v) ~ na.base.i + ka.base.i + input_current
    ]

    return ODESystem(eqs, t, vars, pars; systems=systems, name=name)
end


@component function HHNeuron(;name, input_current=1.0, kwargs...)

    pars = @parameters begin
        Co = 0.01
        Ca = 1
    end

    vars = @variables begin
        v(t) = -70
        tmp_i(t), [input=true]
        #ca_influx() = 0.0
        #i(t), Will be useful and necessary when modeling networks, synapses
    end
    @named constant = Constant(k=input_current)
    @named ground = Ground()
    @named soma = InjectedCapacitor(C=Ca, v=v)
    @named ca_dynamics = CalciumDynamics()

    systems = @named begin
        na = HHSodiumChannel(v_in=v, conductance=100.0, reversal_potential=50.0, parent=ca_dynamics)
        cas = SlowCalciumChannel(v_in=v)
        cat = TransientCalciumChannel(v_in=v)
        ka = ATypePotassiumChannel(v_in=v, conductance=5.0, reversal_potential=-80.0, parent=ca_dynamics)
        kca = CalciumActivatedPotassiumChannel(v_in=v)
        kdr = DelayedRectifierPotassiumChannel(v_in=v, conductance=20.0, reversal_potential=-80.0, parent=ca_dynamics)
        h = HCurrentChannel(v_in=v, conductance=0.5, reversal_potential=-20.0, parent=ca_dynamics)
        leak = LeakChannel(v_in=v, conductance=0.01, reversal_potential=-50.0, parent=ca_dynamics)
    end
    
    systems = [systems..., constant, ground, soma, ca_dynamics]
    eqs = [
        connect(constant.output, soma.I),
        connect(soma.n, ground.g),

        connect(cas.local_ECA, ca_dynamics.E_Ca),
        connect(cat.local_ECA, ca_dynamics.E_Ca),
        connect(kca.local_ECA, ca_dynamics.E_Ca),
        connect(soma.i_fake, tmp_i),
        #connect(v, soma.v),
        #=
        connect(soma.i_output, na.base.i),
        connect(soma.i_output, cas.base.i),
        connect(soma.i_output, cat.base.i),
        connect(soma.i_output, ka.base.i),
        connect(soma.i_output, kca.base.i),
        connect(soma.i_output, kdr.base.i),
        connect(soma.i_output, h.base.i),
        connect(soma.i_output, leak.base.i),=#
        #loop_connect_channels(soma, systems, input_current),
        v ~ soma.v,
        tmp_i ~ summ_channels(systems) + input_current,
        D(ca_dynamics.Ca) ~ (1/ca_dynamics.tau_Ca) * (-ca_dynamics.Ca + ca_dynamics.Ca_base) + 0.94*(cas.base.i + cat.base.i)
    ]
    return ODESystem(eqs, t, vars, pars; systems=systems, name=name)
end

function loop_connect_channels(soma, systems, input_current) #Make into a mltkmodel struct to handle changing input currents
    connections = []
    @variables intensity_sum, [output=true]
    push!(connections, intensity_sum ~ sum(channel.base.i for channel in systems if hasproperty(channel, :base)) + input_current)
    push!(connections, connect(intensity_sum, soma.I))
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