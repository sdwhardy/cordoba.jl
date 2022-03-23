#In DC flow model
function variable_dcgrid_voltage_magnitude_ne(pm::_PM.AbstractDCPModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    # not used
end

#Not used in DC flow model
function variable_dcgrid_voltage_magnitude_ne(pm::_PM.AbstractWModels; kwargs...)
    variable_dcgrid_voltage_magnitude_sqr_ne(pm; kwargs...)
    variable_dcgrid_voltage_magnitude_sqr_du(pm; kwargs...) # duplicated to cancel out existing dc voltages(W) from ohms constraint when z = 0
end

#Not used in DC flow model
function variable_dcgrid_voltage_magnitude_sqr_ne(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    bi_bp = Dict([(i, (b["fbusdc"], b["tbusdc"])) for (i,b) in _PM.ref(pm, nw, :branchdc_ne)])
    bus_vdcmax = merge(Dict([(b,bus["Vdcmax"]) for (b,bus) in _PM.ref(pm, nw, :busdc)]),
    Dict([(b,bus["Vdcmax"]) for (b,bus) in _PM.ref(pm, nw, :busdc_ne)]))
    bus_vdcmin = merge(Dict([(b,bus["Vdcmin"]) for (b,bus) in _PM.ref(pm, nw, :busdc)]),
    Dict([(b,bus["Vdcmin"]) for (b,bus) in _PM.ref(pm, nw, :busdc_ne)]))
         # display(_PM.ids(pm, nw, :buspairsdc_ne))
        wdc_ne = _PM.var(pm, nw)[:wdc_ne] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :busdc)], base_name="$(nw)_wdc_ne",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :busdc, i), "Vdc",  1.0)^2,
        )
        wdcr_ne = _PM.var(pm, nw)[:wdcr_ne] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_wdcr_ne",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :busdc, bi_bp[l][1]), "Vdc",  1.0)^2,
        )
        if bounded
            for (i, busdc) in _PM.ref(pm, nw, :busdc)
                JuMP.set_lower_bound(wdc_ne[i],  busdc["Vdcmin"]^2)
                JuMP.set_upper_bound(wdc_ne[i],  busdc["Vdcmax"]^2)
            end
            for (br, branchdc) in _PM.ref(pm, nw, :branchdc_ne)
                JuMP.set_lower_bound(wdcr_ne[br],  0)
                JuMP.set_upper_bound(wdcr_ne[br],  bus_vdcmax[bi_bp[br][1]] * bus_vdcmax[bi_bp[br][2]])
            end
        end
        report && _IM.sol_component_value(pm, nw, :busdc, :wdc_ne, _PM.ids(pm, nw, :busdc), wdc_ne)
        report && _IM.sol_component_value(pm, nw, :branchdc_ne, :wdcr_ne, _PM.ids(pm, nw, :branchdc_ne), wdcr_ne)
end

function variable_dcgrid_voltage_magnitude_ne(pm::_PM.AbstractLPACModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
        phivdcm_ne = _PM.var(pm, nw)[:phi_vdcm_ne] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :busdc_ne)], base_name="$(nw)_phi_vdcm_ne",
        lower_bound = _PM.ref(pm, nw, :busdc_ne, i, "Vdcmin") - 1,
        upper_bound = _PM.ref(pm, nw, :busdc_ne, i, "Vdcmax") - 1,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :busdc_ne, i), "Vdc")
        )
        if bounded
            for (i, busdc) in _PM.ref(pm, nw, :busdc_ne)
                JuMP.set_lower_bound(phivdcm_ne[i],  busdc["Vdcmin"] - 1)
                JuMP.set_upper_bound(phivdcm_ne[i],  busdc["Vdcmax"] -1 )
            end
        end
        report && _IM.sol_component_value(pm, nw, :busdc_ne, :phivdcm_ne, _PM.ids(pm, nw, :busdc_ne), phivdcm_ne)

#TODO
# think about creating an arc/dict with branchdc_ne pointing to both existing and new buses. Then update limits with corresponding buses.
        phivdcm_fr_ne = _PM.var(pm, nw)[:phi_vdcm_fr] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_phi_vdcm_fr",
        start = 0
        )
        if bounded
            for (i, branchdc) in _PM.ref(pm, nw, :branchdc_ne)
                JuMP.set_lower_bound(phivdcm_fr_ne[i],  -0.2)
                JuMP.set_upper_bound(phivdcm_fr_ne[i],  0.2 )
            end
        end
        report && _IM.sol_component_value(pm, nw, :branchdc_ne, :phivdcm_fr, _PM.ids(pm, nw, :branchdc_ne), phivdcm_fr_ne)


        phivdcm_to_ne = _PM.var(pm, nw)[:phi_vdcm_to] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_phi_vdcm_to",
        start = 0
        )

        if bounded
            for (i, branchdc) in _PM.ref(pm, nw, :branchdc_ne)
                JuMP.set_lower_bound(phivdcm_to_ne[i],  -0.2)
                JuMP.set_upper_bound(phivdcm_to_ne[i],  0.2 )
            end
        end
        report && _IM.sol_component_value(pm, nw, :branchdc_ne, :phivdcm_to, _PM.ids(pm, nw, :branchdc_ne), phivdcm_to_ne)
end

function variable_dcgrid_voltage_magnitude_sqr_du(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true) # this has to to every branch, different than its counterpart(Wdc_fr) since two candidate branches can be connected to same node and two duplicate variables will be needed
    bi_bp = Dict([(i, (b["fbusdc"], b["tbusdc"])) for (i,b) in _PM.ref(pm, nw, :branchdc_ne)])
    wdc_fr_ne = _PM.var(pm, nw)[:wdc_du_fr] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_wdc_du_fr",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :busdc, bi_bp[i][1]), "Vdc",  1.0)^2,
    )
    wdc_to_ne = _PM.var(pm, nw)[:wdc_du_to] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_wdc_du_to",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :busdc, bi_bp[i][1]), "Vdc",  1.0)^2,
    )
    #TODO replace wdc_du_fr and wdc_du_to with wdc_fr and wdc_to make make it consistent with PM, there multiplication is defined by wr - real and wi- imag
    wdcr_frto_ne = _PM.var(pm, nw)[:wdcr_du] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :branchdc_ne)], base_name="$(nw)_wdcr_du",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :busdc, bi_bp[i][1]), "Vdc",  1.0)^2,
    )

    if bounded
        for (i, branchdc) in _PM.ref(pm, nw, :branchdc_ne)
            JuMP.set_lower_bound(wdc_fr_ne[i],  0)
            JuMP.set_upper_bound(wdc_fr_ne[i],  1.21)
            JuMP.set_lower_bound(wdc_to_ne[i],  0)
            JuMP.set_upper_bound(wdc_to_ne[i],  1.21)
            JuMP.set_lower_bound(wdcr_frto_ne[i],  0)
            JuMP.set_upper_bound(wdcr_frto_ne[i],  1.21)
        end
    end
    report && _IM.sol_component_value(pm, nw, :branchdc_ne, :wdc_du_fr, _PM.ids(pm, nw, :branchdc_ne), wdc_fr_ne)
    report && _IM.sol_component_value(pm, nw, :branchdc_ne, :wdc_du_to, _PM.ids(pm, nw, :branchdc_ne), wdc_to_ne)
    report && _IM.sol_component_value(pm, nw, :branchdc_ne, :wdcr_du, _PM.ids(pm, nw, :branchdc_ne), wdcr_frto_ne)
end

#generalized time based constraint on variables (can only expand capacity each year and only once) - for wfs see below
function constraint_t0t1(vss, pm)
    #println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!vdb start: "
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    s=1
    y=1
    for (i,vs) in enumerate(vss)
        for (j,v) in enumerate(last(vs))
            if (mod(i,hl)!=1)
                #println()
                #print(last(vss[i])[j])
                #print("==")
                #print(last(vss[i-1])[j])
                JuMP.@constraint(pm.model, last(vss[i])[j] == last(vss[i-1])[j])
            end
            if (i+hl*yl<=length(vss))
                #println()
                #print(last(vss[i])[j])
                #print("==")
                #print(last(vss[i+hl*yl])[j])
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
end

#time based constraint on wf variables (can only expand capacity each year and only once)
function constraint_t0t1_wfz(vss, pm)
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    s=1
    y=1
    for (i,vs) in enumerate(vss)
        for (j,v) in enumerate(last(vs))
            if (mod(i,hl)!=1)
                JuMP.@constraint(pm.model, last(vss[i])[Int8(j+length(pm.setting["genz"])/2)] == last(vss[i-1])[Int8(j+length(pm.setting["genz"])/2)])
            end
            if (i+hl*yl<=length(vss))
                JuMP.@constraint(pm.model, last(vss[i])[Int8(j+length(pm.setting["genz"])/2)]  == last(vss[i+hl*yl])[Int8(j+length(pm.setting["genz"])/2)])
            end
            if (i==y*hl+(s-1)*yl*hl+1)
                y+=1
            end
            if (i==s*yl*hl+1)
                s+=1
                y=1
            end
            if (i+hl<=s*yl*hl && i+hl<=sl*yl*hl && y<yl)
                    JuMP.@constraint(pm.model, last(vss[i])[Int8(j+length(pm.setting["genz"])/2)]  <= last(vss[i+hl])[Int8(j+length(pm.setting["genz"])/2)])
            end
        end
    end
end

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


#=function constraint_convdc_t0t1(vss, pm)
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
end=#

#=function constraint_branchdc_ne_t0t1(vss, pm)
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
end=#

#=
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
end=#

function variable_dc_converter(pm::_PM.AbstractPowerModel; kwargs...)
    _PMACDC.variable_conv_tranformer_flow(pm; kwargs...)
    _PMACDC.variable_conv_reactor_flow(pm; kwargs...)

    variable_converter_active_power(pm; kwargs...)
    _PMACDC.variable_dcside_power(pm; kwargs...)
    _PMACDC.variable_converter_filter_voltage(pm; kwargs...)
    _PMACDC.variable_converter_internal_voltage(pm; kwargs...)
    _PMACDC.variable_converter_to_grid_active_power(pm; kwargs...)

    _PMACDC.variable_converter_firing_angle(pm; kwargs...)
    _PMACDC.variable_converter_reactive_power(pm; kwargs...)
    _PMACDC.variable_acside_current(pm; kwargs...)
    _PMACDC.variable_converter_to_grid_reactive_power(pm; kwargs...)
end

#Creates required variables for HVDC converter. All are same as PMACDC except active power as it must be less than a variable max  capacity
function variable_dc_converter(pm::_PM.AbstractDCPModel; kwargs...)

    variable_converter_active_power(pm; kwargs...)
    _PMACDC.variable_dcside_power(pm; kwargs...)
    _PMACDC.variable_converter_filter_voltage(pm; kwargs...)
    _PMACDC.variable_converter_internal_voltage(pm; kwargs...)
    _PMACDC.variable_converter_to_grid_active_power(pm; kwargs...)

    _PMACDC.variable_conv_transformer_active_power_to(pm; kwargs...)
    _PMACDC.variable_conv_reactor_active_power_from(pm; kwargs...)
end

#Creates converter max capacity variable and sets  upper/lower limit
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

#constraint on convexified converter AC side power
function variable_converter_active_power(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    pc = _PM.var(pm, nw)[:pconv_ac] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_pconv_ac",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            p_pacmax = _PM.var(pm, nw, :p_pacmax, c)
            JuMP.@constraint(pm.model, pc[c]-p_pacmax  <= 0)
            JuMP.@constraint(pm.model, pc[c]+p_pacmax  >= 0)
        end
    report && _IM.sol_component_value(pm, nw, :convdc, :pconv, _PM.ids(pm, nw, :convdc), pc)
end

###################################### Generators #########################################
#creates variables for WFs
function variable_wfs_peak(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    wf_pacmax = _PM.var(pm, nw)[:wf_pacmax] = JuMP.@variable(pm.model,
    [i in intersect(_PM.ids(pm, nw, :gen),first.(pm.setting["wfz"]))], base_name="$(nw)_wf_pacmax",
    start = 0)

    if bounded
        for (s, gen) in _PM.ref(pm, nw, :gen)
            if issubset([s],first.(pm.setting["wfz"]))
                ########################################
                JuMP.set_lower_bound(wf_pacmax[s],  0)
                JuMP.set_upper_bound(wf_pacmax[s],  last(pm.setting["wfz"][Int8(s-length(pm.setting["genz"])/2)]))
                #######################################
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :wf_pacmax, intersect(_PM.ids(pm, nw, :gen),first.(pm.setting["wfz"])), wf_pacmax)
    return (nw,wf_pacmax)
end

#sets real and imaginary power for a generator
function variable_gen_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_gen_power_real(pm; kwargs...)
    variable_gen_power_imaginary(pm; kwargs...)
end

#real power variable for generators (both conventional and wind)
function variable_gen_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    pg = _PM.var(pm, nw)[:pg] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_pg",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "pg_start")
    )

    if bounded
        for (i, gen) in _PM.ref(pm, nw, :gen)
            if issubset([i],first.(pm.setting["genz"]))
                JuMP.set_lower_bound(pg[i], gen["pmin"])
                JuMP.set_upper_bound(pg[i], gen["pmax"])
            elseif issubset([i],first.(pm.setting["wfz"]))
                wf_pacmax = _PM.var(pm, nw, :wf_pacmax, i)
                JuMP.@constraint(pm.model, pg[i]-gen["pmax"]*wf_pacmax  <= 0)
                JuMP.@constraint(pm.model, pg[i]+gen["pmin"]  >= 0)
            end
        end
    end
    report && _IM.sol_component_value(pm, nw, :gen, :pg, _PM.ids(pm, nw, :gen), pg)
end

#generator imaginary power + constraints
function variable_gen_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    qg = _PM.var(pm, nw)[:qg] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_qg",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "qg_start")
    )

    if bounded
        for (i, gen) in _PM.ref(pm, nw, :gen)
            JuMP.set_lower_bound(qg[i], gen["qmin"])
            JuMP.set_upper_bound(qg[i], gen["qmax"])
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :qg, _PM.ids(pm, nw, :gen), qg)
end

###################################### Storage #########################################
function constraint_storage_losses(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PM.ref(pm, nw, :storage, i)

    _PM.constraint_storage_losses(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
end

function constraint_storage_state(pm::_PM.AbstractPowerModel, i::Int, nw_1::Int, nw_2::Int)
    storage = _PM.ref(pm, nw_2, :storage, i)

    if haskey(_PM.ref(pm, nw_2), :time_elapsed)
        time_elapsed = _PM.ref(pm, nw_2, :time_elapsed)
    else
        Memento.warn(_LOGGER, "network $(nw_2) should specify time_elapsed, using 1.0 as a default")
        time_elapsed = 1.0
    end

    if haskey(_PM.ref(pm, nw_1, :storage), i)
        _FP.constraint_storage_state(pm, nw_1, nw_2, i, storage["charge_efficiency"], storage["discharge_efficiency"], storage["stationary_energy_inflow"], storage["stationary_energy_outflow"], storage["self_discharge_rate"], time_elapsed)
    else
        # if the storage device has status=0 in nw_1, then the stored energy variable will not exist. Initialize storage from data model instead.
        Memento.warn(_LOGGER, "storage component $(i) was not found in network $(nw_1) while building constraint_storage_state between networks $(nw_1) and $(nw_2). Using the energy value from the storage component in network $(nw_2) instead")
        constraint_storage_state_initial(pm, nw_2, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], storage["stationary_energy_inflow"], storage["stationary_energy_outflow"], storage["self_discharge_rate"], time_elapsed)
    end
end

function constraint_storage_thermal_limit(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    ps = _PM.var(pm, nw, :ps, i)
    e_absmax = _PM.var(pm, nw, :e_absmax, i)
    JuMP.@constraint(pm.model, ps-e_absmax/2  <= 0)
    JuMP.@constraint(pm.model, ps+e_absmax/4  >= 0)
end

function constraint_storage_excl_slack(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    sc = _PM.var(pm, nw, :sc, i)
    sd = _PM.var(pm, nw, :sd, i)
    e_absmax = _PM.var(pm, nw, :e_absmax, i)

    JuMP.@constraint(pm.model, sc + sd <= e_absmax/2)
end

function variable_storage_peak(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    e_absmax = _PM.var(pm, nw)[:e_absmax] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_e_absmax",
    start = 0)

    if bounded
        for (s, strg) in _PM.ref(pm, nw, :storage)
            #println("p_rateA[s]: ")
            #println(p_rateA[s])
            if issubset([s],first.(pm.setting["genz"]))
                ########################################
                JuMP.set_lower_bound(e_absmax[s],  0)
                JuMP.set_upper_bound(e_absmax[s],  pm.setting["strg_lim_onshore"])
                #######################################
            elseif issubset([s],first.(pm.setting["wfz"]))
                ########################################
                JuMP.set_lower_bound(e_absmax[s],  0)
                JuMP.set_upper_bound(e_absmax[s],  pm.setting["strg_lim_offshore"])
                #######################################
            end

        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :e_absmax, _PM.ids(pm, nw, :storage), e_absmax)
    #_IM.sol_component_value(pm, nw, :storage, :e_absmax, _PM.ids(pm, nw, :storage), e_absmax)
    return (nw,e_absmax)
end

function constraint_storage_state(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PM.ref(pm, nw, :storage, i)

    if haskey(_PM.ref(pm, nw), :time_elapsed)
        time_elapsed = _PM.ref(pm, nw, :time_elapsed)
    else
        Memento.warn(_LOGGER, "network data should specify time_elapsed, using 1.0 as a default")
        time_elapsed = 1.0
    end
    constraint_storage_state_initial(pm, nw, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], storage["stationary_energy_inflow"], storage["stationary_energy_outflow"], storage["self_discharge_rate"], time_elapsed)
end

function constraint_storage_state_initial(pm::_PM.AbstractPowerModel, n::Int, i::Int, energy, charge_eff, discharge_eff, inflow, outflow, self_discharge_rate, time_elapsed)
    sc = _PM.var(pm, n, :sc, i)
    sd = _PM.var(pm, n, :sd, i)
    se = _PM.var(pm, n, :se, i)
    e_absmax = _PM.var(pm, n, :e_absmax, i)

    JuMP.@constraint(pm.model, se == ((1-self_discharge_rate)^time_elapsed)*e_absmax + time_elapsed*(charge_eff*sc - sd/discharge_eff + inflow - outflow))
end


function variable_storage_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_storage_power_real(pm; kwargs...)
    variable_storage_power_imaginary(pm; kwargs...)
    variable_storage_power_control_imaginary(pm; kwargs...)
    _PM.variable_storage_current(pm; kwargs...)
    variable_storage_energy(pm; kwargs...)
    variable_storage_charge(pm; kwargs...)
    variable_storage_discharge(pm; kwargs...)
end

function variable_storage_discharge(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    sd = _PM.var(pm, nw)[:sd] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_sd",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :storage, i), "sd_start", 1)
    )
    #the discharge rate is full discharge in 2 hours - dr=0.5 in maths
    for (s, strg) in _PM.ref(pm, nw, :storage)
        e_absmax = _PM.var(pm, nw, :e_absmax, s)
        JuMP.@constraint(pm.model, sd[s]-e_absmax/2  <= 0)
        JuMP.@constraint(pm.model, sd[s]  >= 0)
    end

    report && _IM.sol_component_value(pm, nw, :storage, :sd, _PM.ids(pm, nw, :storage), sd)
end

function variable_storage_charge(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    sc = _PM.var(pm, nw)[:sc] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_sc",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :storage, i), "sc_start", 1)
    )
    #the charge rate is full charge in 4 hours - cr=0.25 in maths
    for (s, strg) in _PM.ref(pm, nw, :storage)
        e_absmax = _PM.var(pm, nw, :e_absmax, s)
        JuMP.@constraint(pm.model, sc[s]-e_absmax/4  <= 0)
        JuMP.@constraint(pm.model, sc[s]  >= 0)
    end

    report && _IM.sol_component_value(pm, nw, :storage, :sc, _PM.ids(pm, nw, :storage), sc)
end

function variable_storage_energy(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    se = _PM.var(pm, nw)[:se] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_se",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :storage, i), "se_start", 1)
    )

    for (s, strg) in _PM.ref(pm, nw, :storage)
        e_absmax = _PM.var(pm, nw, :e_absmax, s)
        JuMP.@constraint(pm.model, se[s]-e_absmax  <= 0)
        JuMP.@constraint(pm.model, se[s]  >= 0)
    end

    report && _IM.sol_component_value(pm, nw, :storage, :se, _PM.ids(pm, nw, :storage), se)
end

function variable_storage_power_control_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    qsc = _PM.var(pm, nw)[:qsc] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_qsc",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :storage, i), "qsc_start")
    )

        for (s, strg) in _PM.ref(pm, nw, :storage)
            e_absmax = _PM.var(pm, nw, :e_absmax, s)
            JuMP.@constraint(pm.model, qsc[s]-e_absmax/4  <= 0)
            JuMP.@constraint(pm.model, qsc[s]+e_absmax/4  >= 0)
        end
    report && _IM.sol_component_value(pm, nw, :storage, :qsc, _PM.ids(pm, nw, :storage), qsc)
end

function variable_storage_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    ps = _PM.var(pm, nw)[:ps] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_ps",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :storage, i), "ps_start")
    )

    for (s, strg) in _PM.ref(pm, nw, :storage)
        e_absmax = _PM.var(pm, nw, :e_absmax, s)
        JuMP.@constraint(pm.model, ps[s]-e_absmax/4  <= 0)
        JuMP.@constraint(pm.model, ps[s]+e_absmax/4  >= 0)
    end
    report && _IM.sol_component_value(pm, nw, :storage, :ps, _PM.ids(pm, nw, :storage), ps)
end

function variable_storage_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    qs = _PM.var(pm, nw)[:qs] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_qs",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :storage, i), "qs_start")
    )

    for (s, strg) in _PM.ref(pm, nw, :storage)
        e_absmax = _PM.var(pm, nw, :e_absmax, s)
        JuMP.@constraint(pm.model, qs[s]-e_absmax/4  <= 0)
        JuMP.@constraint(pm.model, qs[s]+e_absmax/4  >= 0)
    end

    report && _IM.sol_component_value(pm, nw, :storage, :qs, _PM.ids(pm, nw, :storage), qs)
end

function variable_absorbed_energy(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    e_abs = _PM.var(pm, nw)[:e_abs] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :storage)], base_name="$(nw)_e_abs",
    start = 0)

    for (s, strg) in _PM.ref(pm, nw, :storage)
        e_absmax = _PM.var(pm, nw, :e_absmax, s)
        JuMP.@constraint(pm.model, e_abs[s]-5000*e_absmax  <= 0)
        JuMP.@constraint(pm.model, e_abs[s]  >= 0)
    end

    report && _IM.sol_component_value(pm, nw, :storage, :e_abs, _PM.ids(pm, nw, :storage), e_abs)
end

############################ Binary DC branch constraints #####################################
"variable: `0 <= branch_ne[l] <= 1` for `l` in `branch`es"
function variable_ne_branch_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_branch_ne = _PM.var(pm, nw)[:branch_ne] = JuMP.@variable(pm.model,
            [l in _PM.ids(pm, nw, :ne_branch)], base_name="$(nw)_branch_ne",
            binary = true,
            start = _PM.comp_start_value(_PM.ref(pm, nw, :ne_branch, l), "branch_tnep_start", 1.0)
        )
    else
        z_branch_ne = _PM.var(pm, nw)[:branch_ne] = JuMP.@variable(pm.model,
            [l in _PM.ids(pm, nw, :ne_branch)], base_name="$(nw)_branch_ne",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = _PM.comp_start_value(_PM.ref(pm, nw, :ne_branch, l), "branch_tnep_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :ne_branch, :built, _PM.ids(pm, nw, :ne_branch), z_branch_ne)
    return (nw,z_branch_ne)
end


# same as that in FlexPlan but variable with nw number is returned for time based constraints
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

############################### constraints for dc branch convexafy ###############################
#Upper limit for convexafied dc branch
function variable_dcbranch_peak(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    p_rateA = _PM.var(pm, nw)[:p_rateA] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :branchdc)], base_name="$(nw)_p_rateA",
    start = 0)
    if bounded
        for (s, branchdc) in _PM.ref(pm, nw, :branchdc)
            #######################################
            JuMP.set_lower_bound(p_rateA[s],  0)
            JuMP.set_upper_bound(p_rateA[s],  pm.setting["ic_lim"])
            #######################################
        end
    end

    report && _IM.sol_component_value(pm, nw, :branchdc, :p_rateA, _PM.ids(pm, nw, :branchdc), p_rateA)
    return (nw,p_rateA)
end

#DC branch flow for continuous transmission line - convex approximation
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
        end
    end
    report && _IM.sol_component_value_edge(pm, nw, :branchdc, :pf, :pt, _PM.ref(pm, nw, :arcs_dcgrid_from), _PM.ref(pm, nw, :arcs_dcgrid_to), p)
end

############################### constraints for ac branch convexafy ###############################
#Upper limit for convexafied dc branch
function variable_acbranch_peak(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool = true, report::Bool=true)
    p_rateAC = _PM.var(pm, nw)[:p_rateAC] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :branch)], base_name="$(nw)_p_rateAC",
    start = 0)
    if bounded
        for (s, branch) in _PM.ref(pm, nw, :branch)
            ########################################
            JuMP.set_lower_bound(p_rateAC[s],  0)
            JuMP.set_upper_bound(p_rateAC[s],  pm.setting["rad_lim"])
            #######################################
        end
    end

    report && _IM.sol_component_value(pm, nw, :branch, :p_rateAC, _PM.ids(pm, nw, :branch), p_rateAC)
    return (nw,p_rateAC)
end
#Upper limit for convexafied ac branch
function variable_branch_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_branch_power_real(pm; kwargs...)
    _PM.variable_branch_power_imaginary(pm; kwargs...)
end


"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_branch_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    p = _PM.var(pm, nw)[:p] = JuMP.@variable(pm.model,
        [(l,i,j) in _PM.ref(pm, nw, :arcs)], base_name="$(nw)_p",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "p_start")
    )

    if bounded
        for arc in _PM.ref(pm, nw, :arcs)
            l,i,j = arc
            p_rateAC = _PM.var(pm, nw, :p_rateAC, l)
            JuMP.@constraint(pm.model, p[arc]-p_rateAC  <= 0)
            JuMP.@constraint(pm.model, p[arc]+p_rateAC  >= 0)
        end
    end

    for (l,branch) in _PM.ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            JuMP.set_start_value(p[f_idx], branch["pf_start"])
        end
        if haskey(branch, "pt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            JuMP.set_start_value(p[t_idx], branch["pt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :pf, :pt, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), p)
end

################################################### Power Balance equation ########################
################################## OBZ #######################
# Constraint template: Power balance constraint including candidate storage
# this is the function 2 of the power balance constraint for a OBZ secanario in cluding storage.
# All nodes can have unique shadow prices
#=function constraint_power_balance_acne_dcne_strg(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = PowerModels.ref(pm, nw, :bus, i)
    bus_arcs = PowerModels.ref(pm, nw, :bus_arcs, i)
    bus_arcs_ne = PowerModels.ref(pm, nw, :ne_bus_arcs, i)
    bus_arcs_dc = PowerModels.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = PowerModels.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = PowerModels.ref(pm, nw, :bus_convs_ac, i)
    bus_convs_ac_ne = PowerModels.ref(pm, nw, :bus_convs_ac_ne, i)
    bus_loads = PowerModels.ref(pm, nw, :bus_loads, i)
    bus_shunts = PowerModels.ref(pm, nw, :bus_shunts, i)
    bus_storage = PowerModels.ref(pm, nw, :bus_storage, i)
    bus_storage_ne = 1#PowerModels.ref(pm, nw, :bus_storage_ne, i)

    pd = Dict(k => PowerModels.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    qd = Dict(k => PowerModels.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    gs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
    constraint_power_balance_acne_dcne_strg(pm, nw, i, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
end=#
function constraint_power_balance_acne_dcne_strg(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    #bus = PowerModels.ref(pm, nw, :bus, i)
    bus_arcs = PowerModels.ref(pm, nw, :bus_arcs, i)
    bus_arcs_ne = PowerModels.ref(pm, nw, :ne_bus_arcs, i)
    bus_arcs_dc = PowerModels.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = PowerModels.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = PowerModels.ref(pm, nw, :bus_convs_ac, i)
    bus_convs_ac_ne = PowerModels.ref(pm, nw, :bus_convs_ac_ne, i)
    bus_loads = PowerModels.ref(pm, nw, :bus_loads, i)
    bus_shunts = PowerModels.ref(pm, nw, :bus_shunts, i)
    bus_storage = PowerModels.ref(pm, nw, :bus_storage, i)
    bus_storage_ne = 1#PowerModels.ref(pm, nw, :bus_storage_ne, i)

    pd = Dict(k => PowerModels.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    qd = Dict(k => PowerModels.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    gs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
    if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
        cost=constraint_power_balance_acne_dcne_strg_admm(pm, nw, i, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
        return cost
    else
        constraint_power_balance_acne_dcne_strg(pm, nw, i, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    end
end
#Power balance constraint including candidate storage
function constraint_power_balance_acne_dcne_strg(pm::_PM.AbstractDCPModel, n::Int, i::Int, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    p = _PM.var(pm, n, :p)
    pg = _PM.var(pm, n, :pg)
    pconv_grid_ac_ne = _PM.var(pm, n, :pconv_tf_fr_ne)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    pconv_ac_ne = _PM.var(pm, n, :pconv_ac_ne)
    p_ne = _PM.var(pm, n, :p_ne)
    ps   = _PM.var(pm, n, :ps)
    #ps_ne   = _PM.var(pm, n, :ps_ne)
    v = 1

    #cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) -sum(ps_ne[s] for s in bus_storage_ne) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    #println(cstr)
    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
    end
end

#=
# NOT complete - Imaginary part needs to be added
function constraint_power_balance_acne_dcne_strg(pm::_PM.AbstractPowerModel, n::Int, i::Int, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    p = _PM.var(pm, n, :p)
    pg = _PM.var(pm, n, :pg)
    pconv_grid_ac_ne = _PM.var(pm, n, :pconv_tf_fr_ne)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    pconv_ac_ne = _PM.var(pm, n, :pconv_ac_ne)
    p_ne = _PM.var(pm, n, :p_ne)
    ps   = _PM.var(pm, n, :ps)
    #ps_ne   = _PM.var(pm, n, :ps_ne)
    v = 1

    #cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) -sum(ps_ne[s] for s in bus_storage_ne) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
    end
end
=#


########################################## Home Market #########################
# Constraint template: Power balance constraint including candidate storage
# this is the function 1 of the power balance constraint for a home market secanario in cluding storage.
# All nodes in the home market must have the same shadow price
#=function constraint_power_balance_acne_dcne_strg_hm(pm::_PM.AbstractPowerModel, is::Set{Int64}; nw::Int=pm.cnw)

    #bus = PowerModels.ref(pm, nw, :bus, i)
    bus_arcs=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs=PowerModels.ref(pm, nw, :bus_arcs, v)
        else;bus_arcs=vcat(bus_arcs,PowerModels.ref(pm, nw, :bus_arcs, v))
        end;end

    bus_arcs_ne=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs_ne=PowerModels.ref(pm, nw, :ne_bus_arcs, v)
        else;bus_arcs_ne=vcat(bus_arcs_ne,PowerModels.ref(pm, nw, :ne_bus_arcs, v))
        end;end

    bus_arcs_dc=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs_dc=PowerModels.ref(pm, nw, :bus_arcs_dc, v)
        else;bus_arcs_dc=vcat(bus_arcs_dc,PowerModels.ref(pm, nw, :bus_arcs_dc, v))
        end;end

    bus_gens=[];for (i,v) in enumerate(is);
        if (i==1);bus_gens=PowerModels.ref(pm, nw, :bus_gens, v)
        else;bus_gens=vcat(bus_gens,PowerModels.ref(pm, nw, :bus_gens, v))
        end;end

    bus_convs_ac=[];for (i,v) in enumerate(is);
        if (i==1);bus_convs_ac=PowerModels.ref(pm, nw, :bus_convs_ac, v)
        else;bus_convs_ac=vcat(bus_convs_ac,PowerModels.ref(pm, nw, :bus_convs_ac, v))
        end;end

    bus_convs_ac_ne=[];for (i,v) in enumerate(is);
        if (i==1);bus_convs_ac_ne=PowerModels.ref(pm, nw, :bus_convs_ac_ne, v)
        else;bus_convs_ac_ne=vcat(bus_convs_ac_ne,PowerModels.ref(pm, nw, :bus_convs_ac_ne, v))
        end;end

    bus_loads=[];for (i,v) in enumerate(is);
        if (i==1);bus_loads=PowerModels.ref(pm, nw, :bus_loads, v)
        else;bus_loads=vcat(bus_loads,PowerModels.ref(pm, nw, :bus_loads, v))
        end;end

    bus_shunts=[];for (i,v) in enumerate(is);
        if (i==1);bus_shunts=PowerModels.ref(pm, nw, :bus_shunts, v)
        else;bus_shunts=vcat(bus_shunts,PowerModels.ref(pm, nw, :bus_shunts, v))
        end;end

    bus_storage=[];for (i,v) in enumerate(is);
        if (i==1);bus_storage=PowerModels.ref(pm, nw, :bus_storage, v)
        else;bus_storage=vcat(bus_storage,PowerModels.ref(pm, nw, :bus_storage, v))
        end;end

    bus_storage_ne=[];#=for (i,v) in enumerate(is);
        if (i==1);bus_storage_ne=PowerModels.ref(pm, nw, :bus_storage_ne, v)
        else;bus_storage_ne=vcat(bus_storage_ne,PowerModels.ref(pm, nw, :bus_storage_ne, v))
        end;end=#
    pd = Dict(k => PowerModels.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    qd = Dict(k => PowerModels.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    gs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
    constraint_power_balance_acne_dcne_strg_hm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
end=#

function constraint_power_balance_acne_dcne_strg_hm(pm::_PM.AbstractPowerModel, is::Set{Int64}; nw::Int=pm.cnw)

    #bus = PowerModels.ref(pm, nw, :bus, i)
    bus_arcs=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs=PowerModels.ref(pm, nw, :bus_arcs, v)
        else;bus_arcs=vcat(bus_arcs,PowerModels.ref(pm, nw, :bus_arcs, v))
        end;end
    #println("bus_arcs")
    #println(bus_arcs)

    bus_arcs_ne=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs_ne=PowerModels.ref(pm, nw, :ne_bus_arcs, v)
        else;bus_arcs_ne=vcat(bus_arcs_ne,PowerModels.ref(pm, nw, :ne_bus_arcs, v))
        end;end
#    println("bus_arcs_ne")
#    println(bus_arcs_ne)

    bus_arcs_dc=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs_dc=PowerModels.ref(pm, nw, :bus_arcs_dc, v)
        else;bus_arcs_dc=vcat(bus_arcs_dc,PowerModels.ref(pm, nw, :bus_arcs_dc, v))
        end;end
#    println("bus_arcs_dc")
#    println(bus_arcs_dc)

    bus_gens=[];for (i,v) in enumerate(is);
        if (i==1);bus_gens=PowerModels.ref(pm, nw, :bus_gens, v)
        else;bus_gens=vcat(bus_gens,PowerModels.ref(pm, nw, :bus_gens, v))
        end;end
#    println("bus_gens")
#    println(bus_gens)

    bus_convs_ac=[];for (i,v) in enumerate(is);
        if (i==1);bus_convs_ac=PowerModels.ref(pm, nw, :bus_convs_ac, v)
        else;bus_convs_ac=vcat(bus_convs_ac,PowerModels.ref(pm, nw, :bus_convs_ac, v))
        end;end
#    println("bus_convs_ac")
#    println(bus_convs_ac)

    bus_convs_ac_ne=[];for (i,v) in enumerate(is);
        if (i==1);bus_convs_ac_ne=PowerModels.ref(pm, nw, :bus_convs_ac_ne, v)
        else;bus_convs_ac_ne=vcat(bus_convs_ac_ne,PowerModels.ref(pm, nw, :bus_convs_ac_ne, v))
        end;end
#    println("bus_convs_ac_ne")
#    println(bus_convs_ac_ne)

    bus_loads=[];for (i,v) in enumerate(is);
        if (i==1);bus_loads=PowerModels.ref(pm, nw, :bus_loads, v)
        else;bus_loads=vcat(bus_loads,PowerModels.ref(pm, nw, :bus_loads, v))
        end;end
#    println("bus_loads")
#    println(bus_loads)

    bus_shunts=[];for (i,v) in enumerate(is);
        if (i==1);bus_shunts=PowerModels.ref(pm, nw, :bus_shunts, v)
        else;bus_shunts=vcat(bus_shunts,PowerModels.ref(pm, nw, :bus_shunts, v))
        end;end
#    println("bus_shunts")
#    println(bus_shunts)

    bus_storage=[];for (i,v) in enumerate(is);
        if (i==1);bus_storage=PowerModels.ref(pm, nw, :bus_storage, v)
        else;bus_storage=vcat(bus_storage,PowerModels.ref(pm, nw, :bus_storage, v))
        end;end
#    println("bus_storage")
#    println(bus_storage)


    bus_storage_ne=[];#=for (i,v) in enumerate(is);
        if (i==1);bus_storage_ne=PowerModels.ref(pm, nw, :bus_storage_ne, v)
        else;bus_storage_ne=vcat(bus_storage_ne,PowerModels.ref(pm, nw, :bus_storage_ne, v))
        end;end=#
    pd = Dict(k => PowerModels.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    qd = Dict(k => PowerModels.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    gs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
    if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
        cost=constraint_power_balance_acne_dcne_strg_hm_admm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
        return cost
    else
        constraint_power_balance_acne_dcne_strg_hm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    end

end

#Power balance constraint including candidate storage
function constraint_power_balance_acne_dcne_strg_hm(pm::_PM.AbstractDCPModel, n::Int, is::Set{Int64}, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    p = _PM.var(pm, n, :p)
    pg = _PM.var(pm, n, :pg)
    pconv_grid_ac_ne = _PM.var(pm, n, :pconv_tf_fr_ne)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    pconv_ac_ne = _PM.var(pm, n, :pconv_ac_ne)
    p_ne = _PM.var(pm, n, :p_ne)
    ps   = _PM.var(pm, n, :ps)
    #ps_ne   = _PM.var(pm, n, :ps_ne)
    v = 1

    #cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) -sum(ps_ne[s] for s in bus_storage_ne) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
#    println(cstr)
    #=if _IM.report_duals(pm)
        for i in is
            _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr
            _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
        end
    end=#
end

function constraint_power_balance_dc_dcne_hm(pm::_PM.AbstractPowerModel, is::Set{Int64}; nw::Int=pm.cnw)

    #bus = PowerModels.ref(pm, nw, :bus, i)
    bus_arcs_dcgrid=[];for (i,v) in enumerate(is);
        if (i==1);bus_arcs_dcgrid=PowerModels.ref(pm, nw, :bus_arcs_dcgrid, v)
        else;bus_arcs=vcat(bus_arcs_dcgrid,PowerModels.ref(pm, nw, :bus_arcs_dcgrid, v))
        end;end
    #=println("bus_arcs_dcgrid")
    for (l,i,j) in bus_arcs_dcgrid; println(l)print(" ")print(i)print(" ")print(j)
    end=#

    bus_arcs_dcgrid_ne=[];for (i,v) in enumerate(is);

        if (i==1);bus_arcs_dcgrid_ne=PowerModels.ref(pm, nw, :bus_arcs_dcgrid_ne, v)
        else;bus_arcs_dcgrid_ne=vcat(bus_arcs_dcgrid_ne,PowerModels.ref(pm, nw, :bus_arcs_dcgrid_ne, v))
        end;end
    #println("bus_arcs_dcgrid_ne")
    bus_arcs_dcgrid_ne=[a for a in bus_arcs_dcgrid_ne if !(issubset([a[2]],is) && issubset([a[3]],is))] #println(string(a[1])*" "*string(a[2])*" "*string(a[3]))
    #end;end
    #println(bus_arcs_dcgrid_ne)

    bus_convs_dc=[];for (i,v) in enumerate(is);
        if (i==1);bus_convs_dc=PowerModels.ref(pm, nw, :bus_convs_dc, v)
        else;bus_convs_dc=vcat(bus_convs_dc,PowerModels.ref(pm, nw, :bus_convs_dc, v))
        end;end
#    println("bus_convs_dc")
#    println(bus_convs_dc)

    bus_convs_dc_ne=[];for (i,v) in enumerate(is);
        if (i==1);bus_convs_dc_ne=PowerModels.ref(pm, nw, :bus_convs_dc_ne, v)
        else;bus_convs_dc_ne=vcat(bus_convs_dc_ne,PowerModels.ref(pm, nw, :bus_convs_dc_ne, v))
        end;end
#    println("bus_convs_dc_ne")
#    println(bus_convs_dc_ne)

    pd=[];for (i,v) in enumerate(is);
        if (i==1);pd=PowerModels.ref(pm, nw, :busdc, v)["Pdc"]
        else;pd=vcat(pd,PowerModels.ref(pm, nw, :busdc, v)["Pdc"])
        end;end
#    println("pd")
#    println(pd)

    if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
        #cost=constraint_power_balance_acne_dcne_strg_hm_admm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
        return cost
    else
        constraint_power_balance_dc_dcne_hm(pm, nw, is, bus_arcs_dcgrid, bus_arcs_dcgrid_ne, bus_convs_dc, bus_convs_dc_ne, pd)
    end
end

#Power balance constraint including candidate storage
function constraint_power_balance_dc_dcne_hm(pm::_PM.AbstractPowerModel, n::Int, is::Set{Int64}, bus_arcs_dcgrid, bus_arcs_dcgrid_ne, bus_convs_dc, bus_convs_dc_ne, pd)
    p_dcgrid = _PM.var(pm, n, :p_dcgrid)
    p_dcgrid_ne = _PM.var(pm, n, :p_dcgrid_ne)
    pconv_dc = _PM.var(pm, n, :pconv_dc)
    pconv_dc_ne = _PM.var(pm, n, :pconv_dc_ne)

    cstr=JuMP.@constraint(pm.model, sum(p_dcgrid[a] for a in bus_arcs_dcgrid) + sum(p_dcgrid_ne[a] for a in bus_arcs_dcgrid_ne) + sum(pconv_dc[c] for c in bus_convs_dc) + sum(pconv_dc_ne[c] for c in bus_convs_dc_ne)  == -1*sum(pd))
#    println(cstr)
end

function constraint_power_balance_dcne_dcne_hm(pm::_PM.AbstractPowerModel, is::Set{Int64}; nw::Int=pm.cnw)

        bus_i=[];for (i,v) in enumerate(is);
            if (i==1);bus_i=PowerModels.ref(pm, nw, :busdc_ne, v)["busdc_i"]
            else;bus_i=vcat(bus_i,PowerModels.ref(pm, nw, :busdc_ne, v)["busdc_i"])
            end;end
#        println("bus_i")
#        println(bus_i)

        bus_arcs_dcgrid_ne=[];for (i,v) in enumerate(bus_i);
            if (i==1);bus_arcs_dcgrid_ne=PowerModels.ref(pm, nw, :bus_arcs_dcgrid_ne, v)
            else;bus_arcs_dcgrid_ne=vcat(bus_arcs_dcgrid_ne,PowerModels.ref(pm, nw, :bus_arcs_dcgrid_ne, v))
            end;end
    #    println("bus_arcs_dcgrid_ne")
    #    println(bus_arcs_dcgrid_ne)

        bus_ne_convs_dc_ne=[];for (i,v) in enumerate(bus_i);
            if (i==1);bus_ne_convs_dc_ne=PowerModels.ref(pm, nw, :bus_ne_convs_dc_ne, v)
            else;bus_ne_convs_dc_ne=vcat(bus_ne_convs_dc_ne,PowerModels.ref(pm, nw, :bus_ne_convs_dc_ne, v))
            end;end
#        println("bus_ne_convs_dc_ne")
#        println(bus_ne_convs_dc_ne)

        pd_ne=[];for (i,v) in enumerate(is);
            if (i==1);pd_ne=PowerModels.ref(pm, nw, :busdc_ne, v)["Pdc"]
            else;pd_ne=vcat(pd_ne,PowerModels.ref(pm, nw, :busdc_ne, v)["Pdc"])
            end;end
#        println("pd_ne")
#        println(pd_ne)

        if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
            #cost=constraint_power_balance_acne_dcne_strg_hm_admm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
            return cost
        else
            if (length(pd_ne)>0)
                constraint_power_balance_dcne_dcne_hm(pm, nw, is, bus_arcs_dcgrid_ne, bus_ne_convs_dc_ne, pd_ne);end
        end
end


function constraint_power_balance_dcne_dcne_hm(pm::_PM.AbstractPowerModel, n::Int, is::Set{Int64}, bus_arcs_dcgrid_ne, bus_ne_convs_dc_ne, pd_ne)
    p_dcgrid_ne = _PM.var(pm, n, :p_dcgrid_ne)
    pconv_dc_ne = _PM.var(pm, n, :pconv_dc_ne)
    xb = _PM.var(pm, n, :branchdc_ne)
    xc = _PM.var(pm, n, :conv_ne)
    cstr = JuMP.@constraint(pm.model, sum(p_dcgrid_ne[a] for a in bus_arcs_dcgrid_ne) + sum(pconv_dc_ne[c] for c in bus_ne_convs_dc_ne)  == -1*sum(pd_ne))
#    println(cstr)
end


##################### node in zone
function constraint_power_balance_dc_dcne_hm_node(pm::_PM.AbstractPowerModel, i, is::Set{Int64}; nw::Int=pm.cnw)

    #bus = PowerModels.ref(pm, nw, :bus, i)
    bus_arcs_dcgrid = PowerModels.ref(pm, nw, :bus_arcs_dcgrid, i)
    bus_arcs_dcgrid_ne = PowerModels.ref(pm, nw, :bus_arcs_dcgrid_ne, i)

    #println("bus_arcs_dcgrid_ne")
    #bus_arcs_dcgrid_ne_inner=sum([pm.setting["balancing_reserve"]/(1-pm.setting["balancing_reserve"])*PowerModels.ref(pm, nw, :branchdc_ne, a[1])["rateA"]*PowerModels.var(pm, nw, :branchdc_ne, a[1]) for a in bus_arcs_dcgrid_ne if (issubset([a[2]],is) && issubset([a[3]],is))]) #println(string(a[1])*" "*string(a[2])*" "*string(a[3]))
    bus_arcs_dcgrid_ne_inner=sum([(1-pm.setting["balancing_reserve"])*PowerModels.ref(pm, nw, :branchdc_ne, a[1])["rateA"]*PowerModels.var(pm, nw, :branchdc_ne, a[1]) for a in bus_arcs_dcgrid_ne if (issubset([a[2]],is) && issubset([a[3]],is))]) #println(string(a[1])*" "*string(a[2])*" "*string(a[3]))

    bus_arcs_dcgrid_ne=[a for a in bus_arcs_dcgrid_ne if !(issubset([a[2]],is) && issubset([a[3]],is))]
    #end;end
    #println("bus_arcs_dcgrid_ne")
    #println(bus_arcs_dcgrid_ne)
    #println("bus_arcs_dcgrid_ne_inner")
    #println(bus_arcs_dcgrid_ne_inner)

    bus_convs_dc = PowerModels.ref(pm, nw, :bus_convs_dc, i)
    bus_convs_dc_ne = PowerModels.ref(pm, nw, :bus_convs_dc_ne, i)
    bus_convs_dc_ne = PowerModels.ref(pm, nw, :bus_convs_dc_ne, i)
    pd=PowerModels.ref(pm, nw, :busdc, i)["Pdc"]

    if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
        #cost=constraint_power_balance_acne_dcne_strg_hm_admm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
        return cost
    else
        constraint_power_balance_dc_dcne_hm_node(pm, nw, is, bus_arcs_dcgrid, bus_arcs_dcgrid_ne, bus_convs_dc, bus_convs_dc_ne, pd, bus_arcs_dcgrid_ne_inner)
    end
end

#Power balance constraint including candidate storage
function constraint_power_balance_dc_dcne_hm_node(pm::_PM.AbstractPowerModel, n::Int, is::Set{Int64}, bus_arcs_dcgrid, bus_arcs_dcgrid_ne, bus_convs_dc, bus_convs_dc_ne, pd, bus_arcs_dcgrid_ne_inner)
    p_dcgrid = _PM.var(pm, n, :p_dcgrid)
    p_dcgrid_ne = _PM.var(pm, n, :p_dcgrid_ne)
    pconv_dc = _PM.var(pm, n, :pconv_dc)
    pconv_dc_ne = _PM.var(pm, n, :pconv_dc_ne)

    cstr1=JuMP.@constraint(pm.model, sum(p_dcgrid[a] for a in bus_arcs_dcgrid) + sum(p_dcgrid_ne[a] for a in bus_arcs_dcgrid_ne) + sum(pconv_dc[c] for c in bus_convs_dc) + sum(pconv_dc_ne[c] for c in bus_convs_dc_ne) - bus_arcs_dcgrid_ne_inner  <= 0)

    cstr2=JuMP.@constraint(pm.model, sum(p_dcgrid[a] for a in bus_arcs_dcgrid) + sum(p_dcgrid_ne[a] for a in bus_arcs_dcgrid_ne) + sum(pconv_dc[c] for c in bus_convs_dc) + sum(pconv_dc_ne[c] for c in bus_convs_dc_ne) + bus_arcs_dcgrid_ne_inner  >= 0)
    #println(cstr1)
    #println(cstr2)
end


function constraint_power_balance_acne_dcne_strg_hm_node(pm::_PM.AbstractPowerModel, i, is::Set{Int64}; nw::Int=pm.cnw)

    #bus = PowerModels.ref(pm, nw, :bus, i)

    bus_arcs = PowerModels.ref(pm, nw, :bus_arcs, i)
    ne_bus_arcs = PowerModels.ref(pm, nw, :ne_bus_arcs, i)
    #println(ne_bus_arcs)
    #ne_bus_arcs_inner=[pm.setting["balancing_reserve"]/(1-pm.setting["balancing_reserve"])*PowerModels.ref(pm, nw, :ne_branch, a[1])["rate_a"]*PowerModels.var(pm, nw, :branch_ne, a[1]) for a in ne_bus_arcs if (issubset([a[2]],is) && issubset([a[3]],is))] #println(string(a[1])*" "*string(a[2])*" "*string(a[3]))
    ne_bus_arcs_inner=[(1-pm.setting["balancing_reserve"])*PowerModels.ref(pm, nw, :ne_branch, a[1])["rate_a"]*PowerModels.var(pm, nw, :branch_ne, a[1]) for a in ne_bus_arcs if (issubset([a[2]],is) && issubset([a[3]],is))] #println(string(a[1])*" "*string(a[2])*" "*string(a[3]))
    if (length(ne_bus_arcs_inner)>0)
        ne_bus_arcs_inner=sum(ne_bus_arcs_inner)
    else
        ne_bus_arcs_inner=0
    end
    #ne_bus_arcs_inner=sum([1 for a in ne_bus_arcs if (issubset([a[2]],is) && issubset([a[3]],is))]) #println(string(a[1])*" "*string(a[2])*" "*string(a[3]))

    ne_bus_arcs=[a for a in ne_bus_arcs if !(issubset([a[2]],is) && issubset([a[3]],is))]
    #println(ne_bus_arcs)
    #println(ne_bus_arcs_inner)
    bus_arcs_dc = PowerModels.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = PowerModels.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = PowerModels.ref(pm, nw, :bus_convs_ac, i)
    bus_convs_ac_ne = PowerModels.ref(pm, nw, :bus_convs_ac_ne, i)
    bus_loads = PowerModels.ref(pm, nw, :bus_loads, i)
    bus_shunts = PowerModels.ref(pm, nw, :bus_shunts, i)
    bus_storage = PowerModels.ref(pm, nw, :bus_storage, i)
    bus_storage_ne=[];
    pd = Dict(k => PowerModels.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    qd = Dict(k => PowerModels.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    gs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => PowerModels.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
    if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
        cost=constraint_power_balance_acne_dcne_strg_hm_admm(pm, nw, is, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
        return cost
    else
        constraint_power_balance_acne_dcne_strg_hm_node(pm, nw, is, bus_arcs, ne_bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs, ne_bus_arcs_inner)
    end

end

#Power balance constraint including candidate storage
function constraint_power_balance_acne_dcne_strg_hm_node(pm::_PM.AbstractDCPModel, n::Int, is::Set{Int64}, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs, ne_bus_arcs_inner)
    p = _PM.var(pm, n, :p)
    pg = _PM.var(pm, n, :pg)
    pconv_grid_ac_ne = _PM.var(pm, n, :pconv_tf_fr_ne)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    pconv_ac_ne = _PM.var(pm, n, :pconv_ac_ne)
    p_ne = _PM.var(pm, n, :p_ne)
    ps   = _PM.var(pm, n, :ps)
    #ps_ne   = _PM.var(pm, n, :ps_ne)
    v = 1

    #cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) -sum(ps_ne[s] for s in bus_storage_ne) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    cstr1=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  - sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) + sum(pd[d] for d in bus_loads) + sum(gs[s] for s in bus_shunts)*v^2 - ne_bus_arcs_inner <= 0)
    cstr2=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  - sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) + sum(pd[d] for d in bus_loads) + sum(gs[s] for s in bus_shunts)*v^2 + ne_bus_arcs_inner >= 0)
#    println(cstr1)
#    println(cstr2)
    #=if _IM.report_duals(pm)
        for i in is
            _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr
            _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
        end
    end=#
end
