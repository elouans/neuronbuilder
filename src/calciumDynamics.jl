@component function CalciumDynamics(;name)
    pars = @parameters begin
        tau_Ca = 20.0 
        Ca_base = 0.05
    end
    
    vars = @variables begin
        Ca(t) = 0.05
        E_Ca(t), [input=true]
    end

    E_Ca ~ (500.0) * (8.6174e-5) * (283.15) * log(max((3000.0 / Ca),0.001))
    eqs = [
        0~0
    ]
    
    return ODESystem(eqs, t, vars, pars; systems=[], name=name)
end