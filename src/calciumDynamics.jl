module calciumDynamics

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: get_eqs
using OrdinaryDiffEq
using ..utils

export CalciumDynamics

@component function CalciumDynamics(;name, v_in)
    @named pin = Pinwithothername()

    pars = @parameters begin
        tau_Ca = 20.0 
        Ca_base = 0.05 
    end
    
    vars = @variables begin
        Ca(t) = 0.05, [connect=Flow, output=true]
        E_Ca(t)
        I_Ca(t) = 0.0, [input=true]
    end
    
    
    eqs = [
        D(Ca) ~ (1/tau_Ca) * (Ca_base + 0.94*(I_Ca) - Ca),
        E_Ca ~ (500.0) * (8.6174e-5) * (283.15) * log((3000.0 / Ca)),
        pin.i ~ I_Ca,
        pin.v ~ E_Ca
    ]
    
    return ODESystem(eqs, t, vars, pars; systems=[pin], name=name)
end

@connector Pinwithothername begin
    v(t)
    i(t), [connect = Flow]
end
end