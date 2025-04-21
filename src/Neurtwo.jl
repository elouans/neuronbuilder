include("calciumDynamics.jl")
include("ionChannels.jl")


@component function NeuronShell(;name, input_current=0.0, kwargs...)

    pars = @parameters begin 
        gap_C = 0.02
    end

    vars = @variables begin
        v(t) = -70
    end

    @named soma = Soma()
end

@component function Soma(;name, kwargs...)
    pars = @parameters begin
        
    end

    vars = @variables begin
        ca_influx() = 0.0
    end

    systems = @named begin
        
    end
end

