#include("utils.jl")

@component function BaseStaticIonChannel(;name, v_in, conductance, reversal_potential, parent=nothing, kwargs...)
    pars = @parameters begin
        g = conductance
        E = reversal_potential
    end


    vars = @variables begin
        i(t) = 0.0, [output=true]
        m(t) = 0.0
        h(t) = 0.0
    end

    eqs = [
        0~0
    ]

    return ODESystem(eqs, t, vars, pars; systems=[], name=name)
end

@component function BaseDynamicIonChannel(;name, v_in, conductance, parent, kwargs...)
    pars = @parameters begin
        g = conductance
    end


    vars = @variables begin
        i(t), [output=true]
        m(t) = 0.0
        h(t) = 0.0
    end

    eqs = [
        0~0
    ]

    return ODESystem(eqs, t, vars, pars; systems=[], name=name)
end

#In theory setting up this way abstracts some of the annoyances of replicating boiler-plate code.
#End objective is to 'private' most variables, exposing only what the subcomponent biologically impacts.
#In theory this is overwritable: extending ODESystem to use a customArgument class exposing either default or everything, 
#Changeable through Global params.

# 0. Sodium Channel
@component function HHSodiumChannel(;name, v_in, conductance=100.0, reversal_potential=50.0, parent, kwargs...)
    @named base = BaseStaticIonChannel(v_in=v_in, conductance=conductance, reversal_potential=reversal_potential, parent=parent)

    # Activation and inactivation functions
    m_inf(v) = 1.0 / (1.0 + exp((v + 25.5) / -5.29))
    h_inf(v) = 1.0 / (1.0 + exp((v + 48.9) / 5.18))
    tau_m(v) = 1.32 - 1.26 / (1 + exp((v + 120.0) / -25.0))
    tau_h(v) = (0.67 / (1.0 + exp((v + 62.9) / -10.0))) * (1.5 + 1.0 / (1.0 + exp((v + 34.9) / 3.6)))
    
    eqs = [
        D(base.m) ~ (m_inf(v_in) - base.m) / tau_m(v_in),
        D(base.h) ~ (h_inf(v_in) - base.h) / tau_h(v_in),
        base.i ~ base.g * base.m^3 * base.h * (base.E - v_in)
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g, base.E]; systems=[base], name=name)
end

# 1. Slow Calcium Channel
@component function SlowCalciumChannel(;name, v_in, conductance=3.0, parent, kwargs...)
    @named base = BaseDynamicIonChannel(v_in=v_in, conductance=conductance, parent=parent)
    
    m_inf(v) = 1.0 / (1.0 + exp((v + 33.0) / -8.1))
    tau_m(v) = 1.4 + 7.0 / (exp((v + 27.0) / 10.0) + exp((v + 70.0) / -13.0))
    
    eqs = [
        D(base.m) ~ (m_inf(v_in) - base.m) / tau_m(v_in),
        base.i ~ base.g * base.m^2 * (parent - v_in),
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g]; systems=[base], name=name)
end

# 2. Transient Calcium Channel
@component function TransientCalciumChannel(;name, v_in, conductance=1.3, parent, kwargs...)
    @named base = BaseDynamicIonChannel(v_in=v_in, conductance=conductance, parent=parent)
    
    m_inf(v) = 1.0 / (1.0 + exp((v + 50.0) / -7.4))
    h_inf(v) = 1.0 / (1.0 + exp((v + 78.0) / 5.0))
    tau_m(v) = 0.44 + 0.15 / (exp((v + 35.0) / 52.0) + exp((v + 35.0) / -50.0))
    tau_h(v) = 22.7 + 0.27 / (exp((v + 55.0) / 7.0) + exp((v + 55.0) / -7.0))

    eqs = [
        D(base.m) ~ (m_inf(v_in) - base.m) / tau_m(v_in),
        D(base.h) ~ (h_inf(v_in) - base.h) / tau_h(v_in),
        base.i ~ base.g * base.m^3 * base.h * (parent - v_in),
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g]; systems=[base], name=name)
end

# 3. A-type Potassium Channel
@component function ATypePotassiumChannel(;name, v_in, conductance=5.0, reversal_potential=-80.0, parent, kwargs...)
    @named base = BaseStaticIonChannel(v_in=v_in, conductance=conductance, reversal_potential=reversal_potential, parent=parent)
    
    m_inf(v) = 1.0 / (1.0 + exp((v + 60.0) / -8.5))
    h_inf(v) = 1.0 / (1.0 + exp((v + 78.0) / 6.0))
    tau_m(v) = 0.185 + 0.5 / (exp((v + 35.8) / 19.7) + exp((v + 79.7) / -12.7))
    tau_h(v) = 0.5 / (exp((v + 66.0) / 11.0) + exp((v + 66.0) / -11.0))
    
    eqs = [
        D(base.m) ~ (m_inf(v_in) - base.m) / tau_m(v_in),
        D(base.h) ~ (h_inf(v_in) - base.h) / tau_h(v_in),
        base.i ~ base.g * base.m^4 * base.h * (base.E - v_in)
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g, base.E]; systems=[base], name=name)
end

# 4. Calcium-activated Potassium Channel
@component function CalciumActivatedPotassiumChannel(;name, v_in, conductance=10.0, reversal_potential=-80.0, parent, kwargs...)

    @named base = BaseStaticIonChannel(v_in=v_in, conductance=conductance, reversal_potential=reversal_potential, parent=parent, )
    
    m_inf(v, ca) = (ca / (ca + 3.0)) / (1.0 + exp((v + 28.3) / -12.6))

    tau_m(v) = 90.3 - 75.1 / (1.0 + exp((v + 46.0) / -22.7))
    
    eqs = [
        D(base.m) ~ (m_inf(v_in, parent) - base.m) / tau_m(v_in),
        base.i ~ base.g * base.m * (base.E - v_in)
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g, base.E]; systems=[base], name=name)
end

# 5. Delayed Rectifier Potassium Channel
@component function DelayedRectifierPotassiumChannel(;name, v_in, conductance=20.0, reversal_potential=-80.0, parent, kwargs...)
    @named base = BaseStaticIonChannel(v_in=v_in, conductance=conductance, reversal_potential=reversal_potential, parent=parent)
    
    m_inf(v) = 1.0 / (1.0 + exp((v + 12.3) / -11.8))
    tau_m(v) = 7.2 - 6.4 / (1.0 + exp((v + 28.3) / -19.2))
    
    eqs = [
        D(base.m) ~ (m_inf(v_in) - base.m) / tau_m(v_in),
        base.i ~ base.g * base.m^4 * (base.E - v_in)
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g, base.E]; systems=[base], name=name)
end

# 6. H-Current Channel
@component function HCurrentChannel(;name, v_in, conductance=0.5, reversal_potential=-20.0, parent, kwargs...)
    @named base = BaseStaticIonChannel(v_in=v_in, conductance=conductance, reversal_potential=reversal_potential, parent=parent)
    
    m_inf(v) = 1.0 / (1.0 + exp((v + 75.0) / 5.5))
    tau_m(v) = 2.0 / (exp((v + 169.7) / 11.6) + exp((v - 26.7) / -14.3))
    
    eqs = [
        D(base.m) ~ (m_inf(v_in) - base.m) / tau_m(v_in),
        base.i ~ base.g * base.m * (base.E - v_in)
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g, base.E]; systems=[base], name=name)
end

# 7. Leak Current
@component function LeakChannel(;name, v_in, conductance=0.01, reversal_potential=-50.0, kwargs...)
    @named base = BaseStaticIonChannel(v_in=v_in, conductance=conductance, reversal_potential=reversal_potential)
    
    eqs = [
        base.i ~ base.g * (base.E - v_in)
    ]
    
    return ODESystem(eqs, t, [base.m, base.h, base.i], [base.g, base.E]; systems=[base], name=name)
end
