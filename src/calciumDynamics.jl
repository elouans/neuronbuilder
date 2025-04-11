
include("utils.jl")

@component function CalciumDynamics(;name, v_in)
    @named pin = Pin()

    pars = @parameters begin
        tau_Ca = 20.0 
        Ca_base = 0.05 
        I_Ca = 0.0
    end
    
    vars = @variables begin
        Ca(t) = 0.05, [input=true]
        E_Ca(t), [output=true]
    end
    
    
    eqs = [
        D(Ca) ~ (1/tau_Ca) * (Ca_base + 0.94*(I_Ca) - Ca),
        E_Ca ~ (500.0) * (8.6174e-5) * (283.15) * log((3000.0 / Ca)),
        pin.v ~ E_Ca,
        pin.i ~ Ca
    ]
    
    return ODESystem(eqs, t, vars, pars; systems=[pin], name=name)
end

