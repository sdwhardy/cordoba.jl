function constraint_dcbranch_t0t1(vss, pm)
    #println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!vdb start: "
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    s=1
    y=1
    for (i,vs) in enumerate(vss)
        for (j,v) in enumerate(last(vs))
            if (mod(i,hl)!=1)
                JuMP.@constraint(pm.model, last(vss[i])[j] == last(vss[i-1])[j])
            end
            if (i+hl*yl<=length(vss))
                JuMP.@constraint(pm.model, last(vss[i])[j]  == last(vss[i+hl*yl])[j])
            end
            if (i==y*hl+(s-1)*yl*hl+1)
                y+=1
            end
            if (i==s*yl*hl+1)
                s+=1
                y=1
            end
            if (i+hl<=s*yl*hl && i+hl<=sl*yl*hl && y<yl)
                    JuMP.@constraint(pm.model, last(vss[i])[j]  <= last(vss[i+hl])[j])
            end
        end
    end

    #println("vdb end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
end

function constraint_convdc_t0t1(vss, pm)
    #println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!vdb start: "
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    s=1
    y=1
    for (i,vs) in enumerate(vss)
        for (j,v) in enumerate(last(vs))
            if (mod(i,hl)!=1)
                JuMP.@constraint(pm.model, last(vss[i])[j] == last(vss[i-1])[j])
            end
            if (i+hl*yl<=length(vss))
                JuMP.@constraint(pm.model, last(vss[i])[j]  == last(vss[i+hl*yl])[j])
            end
            if (i==y*hl+(s-1)*yl*hl+1)
                y+=1
            end
            if (i==s*yl*hl+1)
                s+=1
                y=1
            end
            if (i+hl<=s*yl*hl && i+hl<=sl*yl*hl && y<yl)
                    JuMP.@constraint(pm.model, last(vss[i])[j]  <= last(vss[i+hl])[j])
            end
        end
    end
    #println("vdb end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
end

function constraint_branchdc_ne_t0t1(vss, pm)
    #println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!vdb start: "
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    s=1
    y=1
    for (i,vs) in enumerate(vss)
        for (j,v) in enumerate(last(vs))
            if (mod(i,hl)!=1)
                JuMP.@constraint(pm.model, last(vss[i])[j] == last(vss[i-1])[j])
            end
            if (i+hl*yl<=length(vss))
                JuMP.@constraint(pm.model, last(vss[i])[j]  == last(vss[i+hl*yl])[j])
            end
            if (i==y*hl+(s-1)*yl*hl+1)
                y+=1
            end
            if (i==s*yl*hl+1)
                s+=1
                y=1
            end
            if (i+hl<=s*yl*hl && i+hl<=sl*yl*hl && y<yl)
                    JuMP.@constraint(pm.model, last(vss[i])[j]  <= last(vss[i+hl])[j])
            end
        end
    end
    #println("vdb end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
end

function variable_dc_converter_ne(pm::_PM.AbstractPowerModel; kwargs...)
    _PMACDC.variable_conv_tranformer_flow_ne(pm; kwargs...)
    _PMACDC.variable_conv_reactor_flow_ne(pm; kwargs...)
    _PMACDC.variable_converter_ne(pm; kwargs...)

    _PMACDC.variable_converter_active_power_ne(pm; kwargs...)
    _PMACDC.variable_converter_reactive_power_ne(pm; kwargs...)
    _PMACDC.variable_acside_current_ne(pm; kwargs...)
    _PMACDC.variable_dcside_power_ne(pm; kwargs...)
    # variable_converter_firing_angle_ne(pm; kwargs...)

    _PMACDC.variable_converter_filter_voltage_ne(pm; kwargs...)
    variable_converter_internal_voltage_ne(pm; kwargs...)
    #
    _PMACDC.variable_converter_to_grid_active_power_ne(pm; kwargs...)
    _PMACDC.variable_converter_to_grid_reactive_power_ne(pm; kwargs...)
end

function variable_converter_internal_voltage_ne(pm::_PM.AbstractPowerModel; kwargs...)
    variable_converter_internal_voltage_magnitude_ne(pm; kwargs...)
    _PMACDC.variable_converter_internal_voltage_angle_ne(pm; kwargs...)
end

function variable_converter_internal_voltage_magnitude_ne(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    vmc_ne = _PM.var(pm, nw)[:vmc_ne] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc_ne)], base_name="$(nw)_vmc_ne",
    start = _PM.ref(pm, nw, :convdc_ne, i, "Vtar")
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc_ne)
            JuMP.set_lower_bound(vmc_ne[c], convdc["Vmmin"])
            JuMP.set_upper_bound(vmc_ne[c], convdc["Vmmax"])
        end
    end
    report && _IM.sol_component_value(pm, nw, :convdc_ne, :vmconv, _PM.ids(pm, nw, :convdc_ne), vmc_ne)
end

function variable_dc_converter(pm::_PM.AbstractPowerModel; kwargs...)
    _PMACDC.variable_conv_tranformer_flow(pm; kwargs...)
    _PMACDC.variable_conv_reactor_flow(pm; kwargs...)

    variable_converter_active_power(pm; kwargs...)
    _PMACDC.variable_converter_reactive_power(pm; kwargs...)
    _PMACDC.variable_acside_current(pm; kwargs...)
    _PMACDC.variable_dcside_power(pm; kwargs...)
    _PMACDC.variable_converter_firing_angle(pm; kwargs...)

    _PMACDC.variable_converter_filter_voltage(pm; kwargs...)
    _PMACDC.variable_converter_internal_voltage(pm; kwargs...)

    _PMACDC.variable_converter_to_grid_active_power(pm; kwargs...)
    _PMACDC.variable_converter_to_grid_reactive_power(pm; kwargs...)
end

function variable_dc_converter(pm::_PM.AbstractDCPModel; kwargs...)
    variable_converter_active_power(pm; kwargs...)
    _PMACDC.variable_dcside_power(pm; kwargs...)
    _PMACDC.variable_converter_filter_voltage(pm; kwargs...)
    _PMACDC.variable_converter_internal_voltage(pm; kwargs...)
    _PMACDC.variable_converter_to_grid_active_power(pm; kwargs...)

    _PMACDC.variable_conv_transformer_active_power_to(pm; kwargs...)
    _PMACDC.variable_conv_reactor_active_power_from(pm; kwargs...)
end

function variable_convdc_peak(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    p_pacmax = _PM.var(pm, nw)[:p_pacmax] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_p_pacmax",
    start = 0)

    if bounded
        for (s, convdc) in _PM.ref(pm, nw, :convdc)
            if (nw>1)
                #p_rateA0 = _PM.var(pm, nw-1, :p_rateA, s)
                #println("p_rateA0: ")
                #println(p_rateA0)
            end
            #println("p_rateA[s]: ")
            #println(p_rateA[s])
            ########################################
            JuMP.set_lower_bound(p_pacmax[s],  0)
            JuMP.set_upper_bound(p_pacmax[s],  pm.setting["ic_lim"])
            #######################################
        end
    end

    report && _IM.sol_component_value(pm, nw, :convdc, :p_pacmax, _PM.ids(pm, nw, :convdc), p_pacmax)
    return (nw,p_pacmax)
end

function variable_converter_active_power(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    pc = _PM.var(pm, nw)[:pconv_ac] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_pconv_ac",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )

    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            p_pacmax = _PM.var(pm, nw, :p_pacmax, c)
            JuMP.@constraint(pm.model, pc[c]-p_pacmax  <= 0)
            JuMP.@constraint(pm.model, pc[c]+p_pacmax  >= 0)
        end
    end

    report && _IM.sol_component_value(pm, nw, :convdc, :pconv, _PM.ids(pm, nw, :convdc), pc)
end

function variable_dcbranch_peak(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    p_rateA = _PM.var(pm, nw)[:p_rateA] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :branchdc)], base_name="$(nw)_p_rateA",
    start = 0)

    if bounded
        for (s, branchdc) in _PM.ref(pm, nw, :branchdc)
            if (nw>1)
                #p_rateA0 = _PM.var(pm, nw-1, :p_rateA, s)
                #println("p_rateA0: ")
                #println(p_rateA0)
            end
            #println("p_rateA[s]: ")
            #println(p_rateA[s])
            ########################################
            JuMP.set_lower_bound(p_rateA[s],  0)
            JuMP.set_upper_bound(p_rateA[s],  pm.setting["ic_lim"])
            #######################################
        end
    end

    report && _IM.sol_component_value(pm, nw, :branchdc, :p_rateA, _PM.ids(pm, nw, :branchdc), p_rateA)
    return (nw,p_rateA)
end

function variable_active_dcbranch_flow(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    p = _PM.var(pm, nw)[:p_dcgrid] = JuMP.@variable(pm.model,
    [(l,i,j) in _PM.ref(pm, nw, :arcs_dcgrid)], base_name="$(nw)_pdcgrid",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :branchdc, l), "p_start", 1.0)
    )
    if bounded
        for arc in _PM.ref(pm, nw, :arcs_dcgrid)
            l,i,j = arc
            p_rateA = _PM.var(pm, nw, :p_rateA, l)
            JuMP.@constraint(pm.model, p[arc]-p_rateA  <= 0)
            JuMP.@constraint(pm.model, p[arc]+p_rateA  >= 0)
            #JuMP.set_lower_bound(p[arc], -p_rateA)
            #JuMP.set_upper_bound(p[arc],  p_rateA)
        end
    end
    report && _IM.sol_component_value_edge(pm, nw, :branchdc, :pf, :pt, _PM.ref(pm, nw, :arcs_dcgrid_from), _PM.ref(pm, nw, :arcs_dcgrid_to), p)
end

"variable: `0 <= convdc_ne[c] <= 1` for `c` in `candidate converters"
function variable_branch_ne(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        Z_dc_branch_ne = _PM.var(pm, nw)[:branchdc_ne] = JuMP.@variable(pm.model, #branch_ne is also name in PowerModels, branchdc_ne is candidate branches
        [l in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_branch_ne",
        binary = true,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branchdc_ne, l), "convdc_tnep_start",  0.0)
        )
    else
        Z_dc_branch_ne = _PM.var(pm, nw)[:branchdc_ne] = JuMP.@variable(pm.model, #branch_ne is also name in PowerModels, branchdc_ne is candidate branches
        [l in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_branch_ne",
        lower_bound = 0,
        upper_bound = 1,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branchdc_ne, l), "convdc_tnep_start",  0.0)
        )
    end
    report && _IM.sol_component_value(pm, nw, :branchdc_ne, :isbuilt, _PM.ids(pm, nw, :branchdc_ne), Z_dc_branch_ne)
    return (nw,Z_dc_branch_ne)
end
