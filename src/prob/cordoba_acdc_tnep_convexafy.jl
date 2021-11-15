export cordoba_mp_acdctnepopf_convexafy

""
function cordoba_mp_acdctnepopf_convexafy(data::Dict{String,Any}, model_type::Type, solver; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PM.ref_add_on_off_va_bounds!], setting = s, kwargs...)
    s = setting
    return _PM.run_model(data, model_type, solver, post_cordoba_mp_acdctnepopf_convexafy; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PM.ref_add_on_off_va_bounds!], setting = s, kwargs...)
end

# Here the problem is defined, which is then sent to the solver.
# It is basically a declarion of variables and constraint of the problem

function post_cordoba_mp_acdctnepopf_convexafy(pm::_PM.AbstractPowerModel)
# VARIABLES: defined within PowerModels(ACDC) can directly be used, other variables need to be defined in the according sections of the code: flexible_demand.jl
    vdp=[];vcp=[]
    for n in _PM.nw_ids(pm)
        _PM.variable_bus_voltage(pm; nw = n)
        _PM.variable_gen_power(pm; nw = n)
        _PM.variable_branch_power(pm; nw = n)

        push!(vdp,variable_dcbranch_peak(pm; nw = n))
        variable_active_dcbranch_flow(pm; nw = n)

        push!(vcp,variable_convdc_peak(pm; nw = n))
        variable_dc_converter(pm; nw = n)



        _PMACDC.variable_dcbranch_current(pm; nw = n)
        _PMACDC.variable_dcgrid_voltage_magnitude(pm; nw = n)
    end
    sort!(vdp, by = x -> x[1])
    constraint_dcbranch_t0t1(vdp,pm)
    sort!(vcp, by = x -> x[1])
    constraint_convdc_t0t1(vcp,pm)
#OBJECTIVE see objective.jl
    objective_min_cost_acdc_convexafy(pm)
#CONSTRAINTS: defined within PowerModels(ACDC) can directly be used, other constraints need to be defined in the according sections of the code: flexible_demand.jl
    for n in _PM.nw_ids(pm)
        _PM.constraint_model_voltage(pm; nw = n)
        _PM.constraint_ne_model_voltage(pm; nw = n)
        _PMACDC.constraint_voltage_dc(pm; nw = n)
        _PMACDC.constraint_voltage_dc_ne(pm; nw = n)
        for i in _PM.ids(pm, n, :ref_buses)
            _PM.constraint_theta_ref(pm, i, nw = n)
        end

        for i in _PM.ids(pm, n, :bus)
            _PMACDC.constraint_power_balance_ac(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :branch)
            _PM.constraint_ohms_yt_from(pm, i; nw = n)
            _PM.constraint_ohms_yt_to(pm, i; nw = n)
            _PM.constraint_voltage_angle_difference(pm, i; nw = n)
            _PM.constraint_thermal_limit_from(pm, i; nw = n)
            _PM.constraint_thermal_limit_to(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :busdc)
            _PMACDC.constraint_power_balance_dc(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :branchdc)
            _PMACDC.constraint_ohms_dc_branch(pm, i; nw = n)
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
    end
    #println(pm.model)
end
