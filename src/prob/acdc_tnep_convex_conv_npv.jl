export acdc_tnep_convex_conv_npv

""
function acdc_tnep_convex_conv_npv(data::Dict{String,Any}, model_type::Type, solver; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.add_candidate_dcgrid!, _PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!], setting = s, kwargs...)
    s = setting
    return _PM.run_model(data, model_type, solver, post_acdc_tnep_convex_conv_npv; ref_extensions = [_PMACDC.add_ref_dcgrid!,_PMACDC.add_candidate_dcgrid!, _PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!], setting = s, kwargs...)
end

# Here the problem is defined, which is then sent to the solver.
# It is basically a declarion of variables and constraint of the problem

""
function post_acdc_tnep_convex_conv_npv(pm::_PM.AbstractPowerModel)
    # VARIABLES: defined within PowerModels(ACDC) can directly be used, other variables need to be defined in the according sections of the code: flexible_demand.jl
    vcp=[];vbp=[]
    for n in _PM.nw_ids(pm)
        _PM.variable_bus_voltage(pm; nw = n)#CREATES voltage angle variables and sets voltage magnitudes to 1.0 at all AC nodes
        _PM.variable_gen_power(pm; nw = n)#CREATES pg variables of all generators and sets upper and lower limits
        _PM.variable_branch_power(pm; nw = n)#CREATES real and imaginary branch power variables, sets upper and lower limits and initial values

        #_PMACDC.variable_dc_converter(pm; nw = n)#CREATES pconv_ac, pconv_dc, vaf, vac(internal converter angle), [pconv_tf_fr, pconv_tf_to(transformer)]-AC side, pconv_pr_fr (phase reactor),
        ####################### cordoba ########################
        push!(vcp,variable_convdc_peak(pm; nw = n))
        variable_dc_converter(pm; nw = n)
        ########################################################

        ######################### May use in future ###########################################
        #_PM.variable_storage_power(pm; nw = n)
        _PMACDC.variable_voltage_slack(pm; nw = n)#CREATES va_du(angle @ transformer?), vaf_du(angle @ filter), vac_du(angle @ converter?) for each convdc_ne and bounds between -2pi and 2pi, intiallizes to 0
        _PMACDC.variable_active_dcbranch_flow(pm; nw = n)#CREATES pdcgrid for each branchdc and bounds between +/- rateA

        _PMACDC.variable_dcbranch_current(pm; nw = n)#ONLY ACTIVE IN BF MODEL(not doing anything even with branch flow = true?), CREATES: ccm_dcgrid and bounds.
        _PMACDC.variable_dcgrid_voltage_magnitude(pm; nw = n)# IN DC flow MODEL not USED
        #variable_absorbed_energy(pm; nw = n)#storage (not investigated yet)
        #variable_absorbed_energy_ne(pm; nw = n)#storage (not investigated yet)
        #variable_flexible_demand(pm; nw = n)#storage (not investigated yet)
        #######################################################################################

        # variables for TNEP problem

        _PM.variable_ne_branch_indicator(pm; nw = n)#CREATES: branch_ne {0,1} for each branch_ne
        _PM.variable_ne_branch_power(pm; nw = n)#CREATES p_ne branch power variable and sets bounds
        push!(vbp,variable_branch_ne(pm; nw = n))#CREATES: branch_ne {0,1} for each branchdc_ne
        _PMACDC.variable_active_dcbranch_flow_ne(pm; nw = n)#CREATES: pdcgrid_ne for each branchdc_ne and sets upper and lower bounds

        ######################### May use in future ###########################################
        _PM.variable_ne_branch_voltage(pm; nw = n)#DOES nothing in DC model - there are no voltages on branches
        #variable_storage_power_ne(pm; nw = n)#storage (not investigated yet)
        _PMACDC.variable_dc_converter_ne(pm; nw = n)#CREATES: pconv_ac_ne, pconv_dc_ne, vaf_ne, vac_ne, pconv_tf_fr_ne, pconv_tf_to_ne, pconv_pr_from_ne and sets bounds (equivalent to existing conv vars)
        _PMACDC.variable_dcbranch_current_ne(pm; nw = n)##ONLY ACTIVE IN BF MODEL(not doing anything even with branch flow = true?),CREATES: ccm_dcgrid_ne
        variable_dcgrid_voltage_magnitude_ne(pm; nw = n)# IN DC flow MODEL not USED=#
        #######################################################################################
    end
    sort!(vbp, by = x -> x[1])
    constraint_t0t1(vbp,pm)
    sort!(vcp, by = x -> x[1])
    constraint_t0t1(vcp,pm)

#CONSTRAINTS: defined within PowerModels(ACDC) can directly be used, other constraints need to be defined in the according sections of the code: flexible_demand.jl
    for n in _PM.nw_ids(pm)
        _PM.constraint_model_voltage(pm; nw = n)
        _PM.constraint_ne_model_voltage(pm; nw = n)
        _PMACDC.constraint_voltage_dc(pm; nw = n)
        _PMACDC.constraint_voltage_dc_ne(pm; nw = n)
        for i in _PM.ids(pm, n, :ref_buses)
            _PM.constraint_theta_ref(pm, i, nw = n)
        end

        #=for i in _PM.ids(pm, n, :bus)
            constraint_power_balance_acne_dcne_flex(pm, i; nw = n)
        end=#
        for i in _PM.ids(pm, n, :bus)
            _PMACDC.constraint_power_balance_acne_dcne(pm, i; nw = n)
        end

        if haskey(pm.setting, "allow_line_replacement") && pm.setting["allow_line_replacement"] == true
            for i in _PM.ids(pm, n, :branch)
                constraint_ohms_yt_from_repl(pm, i; nw = n)
                constraint_ohms_yt_to_repl(pm, i; nw = n)
                constraint_voltage_angle_difference_repl(pm, i; nw = n)
                constraint_thermal_limit_from_repl(pm, i; nw = n)
                constraint_thermal_limit_to_repl(pm, i; nw = n)
            end
        else
            for i in _PM.ids(pm, n, :branch)
                _PM.constraint_ohms_yt_from(pm, i; nw = n)
                _PM.constraint_ohms_yt_to(pm, i; nw = n)
                _PM.constraint_voltage_angle_difference(pm, i; nw = n)
                _PM.constraint_thermal_limit_from(pm, i; nw = n)
                _PM.constraint_thermal_limit_to(pm, i; nw = n)
            end
        end
        for i in _PM.ids(pm, n, :ne_branch)
            _PM.constraint_ne_ohms_yt_from(pm, i; nw = n)
            _PM.constraint_ne_ohms_yt_to(pm, i; nw = n)
            _PM.constraint_ne_voltage_angle_difference(pm, i; nw = n)
            _PM.constraint_ne_thermal_limit_from(pm, i; nw = n)
            _PM.constraint_ne_thermal_limit_to(pm, i; nw = n)
            if n > 1
                _PMACDC.constraint_candidate_acbranches_mp(pm, n, i)
            end
        end

        for i in _PM.ids(pm, n, :busdc)
            _PMACDC.constraint_power_balance_dc_dcne(pm, i; nw = n)
        end
        for i in _PM.ids(pm, n, :busdc_ne)
            _PMACDC.constraint_power_balance_dcne_dcne(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :branchdc)
            _PMACDC.constraint_ohms_dc_branch(pm, i; nw = n)
        end

        candidates_in_corridor=[]
        for i in _PM.ids(pm, n, :branchdc_ne)
            _PMACDC.constraint_ohms_dc_branch_ne(pm, i; nw = n)
            _PMACDC.constraint_branch_limit_on_off(pm, i; nw = n)
            if n > 1
                _PMACDC.constraint_candidate_dcbranches_mp(pm, n, i)
            end
            push!(candidates_in_corridor,collect_4_constraint_candidate_corridor_limit(pm, i; nw = n))
        end

        if (pm.setting["corridor_limit"])
            constraint_candidate_corridor_limit(pm, candidates_in_corridor; nw = n)
        end

        for i in _PM.ids(pm, n, :convdc)
            _PMACDC.constraint_converter_losses(pm, i; nw = n)
            _PMACDC.constraint_converter_current(pm, i; nw = n)
            _PMACDC.constraint_conv_transformer(pm, i; nw = n)
            _PMACDC.constraint_conv_reactor(pm, i; nw = n)
            _PMACDC.constraint_conv_filter(pm, i; nw = n)
            if pm.ref[:nw][n][:convdc][i]["islcc"] == 1
                _PMACDC.constraint_conv_firing_angle(pm, i; nw = n)
            end
        end

        for i in _PM.ids(pm, n, :convdc_ne)
            _PMACDC.constraint_converter_losses_ne(pm, i; nw = n)
            _PMACDC.constraint_converter_current_ne(pm, i; nw = n)
            _PMACDC.constraint_converter_limit_on_off(pm, i; nw = n)
            if n > 1
                _PMACDC.constraint_candidate_converters_mp(pm, n, i)
            end
            _PMACDC.constraint_conv_transformer_ne(pm, i; nw = n)
            _PMACDC.constraint_conv_reactor_ne(pm, i; nw = n)
            _PMACDC.constraint_conv_filter_ne(pm, i; nw = n)
            if pm.ref[:nw][n][:convdc_ne][i]["islcc"] == 1
                _PMACDC.constraint_conv_firing_angle_ne(pm, i; nw = n)
            end
        end
        if (haskey(pm.setting,"agent") && pm.setting["agent"]=="fixed_cables_cons")
            for n in _PM.nw_ids(pm)
                fix_variables(pm, n)
            end
        end
    end
    #OBJECTIVE see objective.jl
    objective_min_cost_acdc_convex_conv_npv(pm)
    #println(JuMP.objective_function(pm.model))
end
