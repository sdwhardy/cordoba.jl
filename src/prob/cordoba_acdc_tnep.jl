export cordoba_mp_acdctnepopf

""
function cordoba_mp_acdctnepopf(file::String, model_type::Type, solver; kwargs...)
    data = _PM.parse_file(file)
    return cordoba_mp_acdctnepopf(data, model_type, solver; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.add_candidate_dcgrid!], kwargs...)
end

""
function cordoba_mp_acdctnepopf(data::Dict{String,Any}, model_type::Type, solver; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.add_candidate_dcgrid!, _PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!], setting = s, kwargs...)
    if setting["process_data_internally"] == true
        # PowerModelsACDC.process_additional_data!(data)
        _PMACDC.process_additional_data!(data)
    end
    s = setting
    return _PM.run_model(data, model_type, solver, post_cordoba_mp_acdctnepopf; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.add_candidate_dcgrid!, _PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!], setting = s, kwargs...)
    # pm = _PM.build_model(data, model_type, post_mp_tnepopf; setting = s, kwargs...)
    # return _PM.optimize_model!(pm, solver; solution_builder = get_solution_acdc_ne)
end

""
function post_cordoba_mp_acdctnepopf(pm::_PM.AbstractPowerModel)
    for (n, networks) in pm.ref[:nw]
        _PM.variable_bus_voltage(pm; nw = n)
        _PM.variable_gen_power(pm; nw = n)
        _PM.variable_branch_power(pm; nw = n)

        _PMACDC.variable_voltage_slack(pm; nw = n)
        _PMACDC.variable_active_dcbranch_flow(pm; nw = n)
        _PMACDC.variable_dc_converter(pm; nw = n)
        _PMACDC.variable_dcbranch_current(pm; nw = n)
        _PMACDC.variable_dcgrid_voltage_magnitude(pm; nw = n)

        # new variables for TNEP problem
        _PM.variable_ne_branch_indicator(pm; nw = n)
        _PM.variable_ne_branch_power(pm; nw = n)
        _PM.variable_ne_branch_voltage(pm; nw = n)
        _PMACDC.variable_active_dcbranch_flow_ne(pm; nw = n)
        _PMACDC.variable_branch_ne(pm; nw = n)
        _PMACDC.variable_dc_converter_ne(pm; nw = n)
        _PMACDC.variable_dcbranch_current_ne(pm; nw = n)
        _PMACDC.variable_dcgrid_voltage_magnitude_ne(pm; nw = n)
    end
    _PMACDC.objective_min_cost_acdc(pm)
    for (n, networks) in pm.ref[:nw]
        _PM.constraint_model_voltage(pm; nw = n)
        _PM.constraint_ne_model_voltage(pm; nw = n)
        _PMACDC.constraint_voltage_dc(pm; nw = n)
        _PMACDC.constraint_voltage_dc_ne(pm; nw = n)
        for i in _PM.ids(pm, n, :ref_buses)
            _PM.constraint_theta_ref(pm, i, nw = n)
        end

        for i in _PM.ids(pm, n, :bus)
            _PMACDC.constraint_power_balance_acne_dcne(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :branch)
            _PM.constraint_ohms_yt_from(pm, i; nw = n)
            _PM.constraint_ohms_yt_to(pm, i; nw = n)
            _PM.constraint_voltage_angle_difference(pm, i; nw = n)
            _PM.constraint_thermal_limit_from(pm, i; nw = n)
            _PM.constraint_thermal_limit_to(pm, i; nw = n)
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
        for i in _PM.ids(pm, :branchdc_ne)
            _PMACDC.constraint_ohms_dc_branch_ne(pm, i; nw = n)
            _PMACDC.constraint_branch_limit_on_off(pm, i; nw = n)
            if n > 1
                _PMACDC.constraint_candidate_dcbranches_mp(pm, n, i)
            end
        end

        for i in _PM.ids(pm, :convdc)
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
            #display(pm.ref)
            if pm.ref[:nw][n][:convdc_ne][i]["islcc"] == 1
                _PMACDC.constraint_conv_firing_angle_ne(pm, i; nw = n)
            end
        end
    end
end
