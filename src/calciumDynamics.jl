@component function CalciumDynamics(;name)
    pars = @parameters begin
        tau_Ca = 20.0 
        Ca_base = 0.05
        I_Ca = 0.0
    end
    
    vars = @variables begin
        Ca(t) = Ca_base
    end
    
    # Create an observed variable for E_Ca
    @variables E_Ca(t)
    
    eqs = [
        D(Ca) ~ (1/tau_Ca) * (Ca_base + 0.94*(I_Ca) - Ca)
    ]
    
    observed = [
        E_Ca ~ (500.0) * (8.6174e-5) * (283.15) * log((3000.0 / Ca))
    ]
    
    return ODESystem(eqs, t, vars, pars; observed=observed, systems=[], name=name)
end