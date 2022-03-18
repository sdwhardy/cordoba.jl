############################# cordoba_acdc_wf_strg_admm ###############
export cordoba_acdc_wf_strg_admm

""
function cordoba_acdc_wf_strg_admm(data::Dict{String,Any}, model_type::Type, solver; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.add_candidate_dcgrid!, _FP.add_candidate_storage!, _PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!], setting = s, kwargs...)
    s = setting
    return _PM.run_model(data, model_type, solver, post_cordoba_acdc_wf_strg_admm; ref_extensions = [_PMACDC.add_ref_dcgrid!,_PMACDC.add_candidate_dcgrid!, _FP.add_candidate_storage!, _PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!], setting = s, kwargs...)
end

# Here the problem is defined, which is then sent to the solver.
# It is basically a declarion of variables and constraint of the problem



""
function post_cordoba_acdc_wf_strg_admm(pm::_PM.AbstractPowerModel)
    # VARIABLES: defined within PowerModels(ACDC) can directly be used, other variables need to be defined in the according sections of the code: flexible_demand.jl
    vsp=[];vcp=[];vap=[];vbp=[];vgp=[];vdp=[];vacp=[]
        for n in _PM.nw_ids(pm)


            _PM.variable_bus_voltage(pm; nw = n)#CREATES voltage angle variables and sets voltage magnitudes to 1.0 at all AC nodes
            push!(vgp,variable_wfs_peak(pm; nw = n))
            variable_gen_power(pm; nw = n)#CREATES pg variables of all generators and sets upper and lower limits

            ####################### cordoba converters ########################
            push!(vcp,variable_convdc_peak(pm; nw = n))

            variable_dc_converter(pm; nw = n)

            ####################### cordoba storage ########################
            push!(vsp,variable_storage_peak(pm; nw = n))
            variable_storage_power(pm; nw = n)#CREATES variables ps, qsc, se, sc, sd and sets upper/lower bounds
            variable_absorbed_energy(pm; nw = n)#e_abs upper/lower limits
            #_FP.variable_absorbed_energy_ne(pm; nw = n)
            #_FP.variable_storage_power_ne(pm; nw = n)
            ########################################################

            ######################### May use in future ###########################################
            if (haskey(pm.setting,"relax_problem") && pm.setting["relax_problem"]==true)
                push!(vacp,variable_acbranch_peak(pm; nw = n))
                variable_branch_power(pm; nw = n)#CREATES real and imaginary branch power variables, sets upper and lower limits and initial values
                push!(vdp,variable_dcbranch_peak(pm; nw = n))
                variable_active_dcbranch_flow(pm; nw = n)
            else
                _PM.variable_branch_power(pm; nw = n)#CREATES real and imaginary branch power variables, sets upper and lower limits and initial values
                _PMACDC.variable_active_dcbranch_flow(pm; nw = n)#CREATES pdcgrid for each branchdc and bounds between +/- rateA
            end
            _PMACDC.variable_voltage_slack(pm; nw = n)#CREATES va_du(angle @ transformer?), vaf_du(angle @ filter), vac_du(angle @ converter?) for each convdc_ne and bounds between -2pi and 2pi, intiallizes to 0
            _PMACDC.variable_dcbranch_current(pm; nw = n)#ONLY ACTIVE IN BF MODEL(not doing anything even with branch flow = true?), CREATES: ccm_dcgrid and bounds.
            _PMACDC.variable_dcgrid_voltage_magnitude(pm; nw = n)# IN DC flow MODEL not USED

            #######################################################################################

            # variables for TNEP problem

            push!(vap,variable_ne_branch_indicator(pm; nw = n))#CREATES: branch_ne {0,1} for each branch_ne
            _PM.variable_ne_branch_power(pm; nw = n)#CREATES p_ne branch power variable and sets bounds
            push!(vbp,variable_branch_ne(pm; nw = n))#CREATES: branch_ne {0,1} for each branchdc_ne
            _PMACDC.variable_active_dcbranch_flow_ne(pm; nw = n)#CREATES: pdcgrid_ne for each branchdc_ne and sets upper and lower bounds

            ######################### May use in future ###########################################

            _PM.variable_ne_branch_voltage(pm; nw = n)#DOES nothing in DC model - there are no voltages on branches
            _PMACDC.variable_dc_converter_ne(pm; nw = n)#CREATES: pconv_ac_ne, pconv_dc_ne, vaf_ne, vac_ne, pconv_tf_fr_ne, pconv_tf_to_ne, pconv_pr_from_ne and sets bounds (equivalent to existing conv vars)
            _PMACDC.variable_dcbranch_current_ne(pm; nw = n)##ONLY ACTIVE IN BF MODEL(not doing anything even with branch flow = true?),CREATES: ccm_dcgrid_ne
            variable_dcgrid_voltage_magnitude_ne(pm; nw = n)# IN DC flow MODEL not USED

            #######################################################################################
        end
        if (haskey(pm.setting,"relax_problem") && pm.setting["relax_problem"]==true)
            sort!(vdp, by = x -> x[1])
            constraint_t0t1(vdp,pm)
            sort!(vacp, by = x -> x[1])
            constraint_t0t1(vacp,pm)
        end

        sort!(vsp, by = x -> x[1])
        constraint_t0t1(vsp,pm)
        sort!(vbp, by = x -> x[1])
        constraint_t0t1(vbp,pm)
        sort!(vap, by = x -> x[1])
        constraint_t0t1(vap,pm)
        sort!(vcp, by = x -> x[1])
        constraint_t0t1(vcp,pm)
        sort!(vgp, by = x -> x[1])
        constraint_t0t1_wfz(vgp,pm)

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
            if (!(haskey(pm.setting, "agent")) || (pm.setting["agent"] == ""))
                for i in _PM.ids(pm, n, :bus)
                    if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
                        if !(issubset(i,pm.setting["home_market"]))
                            constraint_power_balance_acne_dcne_strg(pm, i; nw = n)
                        else
                            is=intersect(_PM.ids(pm, n, :busdc),first.(pm.setting["home_market"]))
                            constraint_power_balance_acne_dcne_strg_hm_node(pm, i, is; nw = n)
                        end
                    else
                        constraint_power_balance_acne_dcne_strg(pm, i; nw = n)
                    end
                end
            end
            if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
                is=intersect(_PM.ids(pm, n, :bus),first.(pm.setting["home_market"]))
                constraint_power_balance_acne_dcne_strg_hm(pm, is; nw = n)
            end

            for i in _PM.ids(pm, n, :branch)
                _PM.constraint_ohms_yt_from(pm, i; nw = n)
                _PM.constraint_ohms_yt_to(pm, i; nw = n)
                constraint_ohms_ac_branch(pm, i; nw = n)#NOTE change for Zonal market
                _PM.constraint_voltage_angle_difference(pm, i; nw = n)
                _PM.constraint_thermal_limit_from(pm, i; nw = n)
                _PM.constraint_thermal_limit_to(pm, i; nw = n)
            end

            ac_candidates_in_corridor=[]
            for i in _PM.ids(pm, n, :ne_branch)
                _PM.constraint_ne_ohms_yt_from(pm, i; nw = n)
                _PM.constraint_ne_ohms_yt_to(pm, i; nw = n)
                constraint_ohms_ac_branch_ne(pm, i; nw = n)#NOTE change for Zonal market
                _PM.constraint_ne_voltage_angle_difference(pm, i; nw = n)
                _PM.constraint_ne_thermal_limit_from(pm, i; nw = n)
                _PM.constraint_ne_thermal_limit_to(pm, i; nw = n)
                push!(ac_candidates_in_corridor,collect_4_constraint_candidate_corridor_limit_ac(pm, i; nw = n))
            end

            for i in _PM.ids(pm, n, :busdc)
                if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
                    if !(issubset(i,pm.setting["home_market"]))
                        _PMACDC.constraint_power_balance_dc_dcne(pm, i; nw = n)
                    else
                        is=intersect(_PM.ids(pm, n, :busdc),first.(pm.setting["home_market"]))
                        constraint_power_balance_dc_dcne_hm_node(pm, i, is; nw = n)
                    end
                else
                    _PMACDC.constraint_power_balance_dc_dcne(pm, i; nw = n)
                end
            end

            if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
                is=intersect(_PM.ids(pm, n, :busdc),first.(pm.setting["home_market"]))
                constraint_power_balance_dc_dcne_hm(pm, is; nw = n)
            end

            for i in _PM.ids(pm, n, :busdc_ne)
                if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
                    if !(issubset(i,pm.setting["home_market"]))
                        _PMACDC.constraint_power_balance_dcne_dcne(pm, i; nw = n)
                    end
                else
                    _PMACDC.constraint_power_balance_dcne_dcne(pm, i; nw = n)
                end
            end

            if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
                is=intersect(_PM.ids(pm, n, :busdc_ne),first.(pm.setting["home_market"]))
                constraint_power_balance_dcne_dcne_hm(pm, is; nw = n)
            end

            for i in _PM.ids(pm, n, :branchdc)
                constraint_ohms_dc_branch(pm, i; nw = n)#NOTE change for Zonal market
            end

            dc_candidates_in_corridor=[]
            for i in _PM.ids(pm, n, :branchdc_ne)
                constraint_ohms_dc_branch_ne(pm, i; nw = n)#NOTE change for Zonal market
                _PMACDC.constraint_branch_limit_on_off(pm, i; nw = n)
                push!(dc_candidates_in_corridor,collect_4_constraint_candidate_corridor_limit_dc(pm, i; nw = n))
            end

            if (pm.setting["corridor_limit"])
                if (length(ac_candidates_in_corridor)>0)
                    constraint_candidate_corridor_limit(pm, ac_candidates_in_corridor; nw = n);end
                if (length(dc_candidates_in_corridor)>0)
                constraint_candidate_corridor_limit(pm, dc_candidates_in_corridor; nw = n);end
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

            for i in _PM.ids(pm, :storage, nw=n)
                constraint_storage_excl_slack(pm, i, nw = n)
                constraint_storage_thermal_limit(pm, i, nw = n)
                constraint_storage_losses(pm, i, nw = n)
            end
        end

        for (s, scenario) in pm.ref[:scenario]

            network_ids = sort(collect(n for (sc, n) in scenario))
            n_1 = network_ids[1]
            n_last = network_ids[end]

            # NW = 1
            for i in _PM.ids(pm, :storage, nw = n_1)
                _FP.constraint_storage_state(pm, i, nw = n_1)
                _FP.constraint_maximum_absorption(pm, i, nw = n_1)
            end

            # NW = last
            for i in _PM.ids(pm, :storage, nw = n_last)
                _FP.constraint_storage_state_final(pm, i, nw = n_last)
            end

            for n_2 in network_ids[2:end]
                for i in _PM.ids(pm, :storage, nw = n_2)
                    _FP.constraint_storage_state(pm, i, n_1, n_2)
                    _FP.constraint_maximum_absorption(pm, i, n_1, n_2)
                end
                n_1 = n_2
            end
        end
        max_investment_per_year(pm)
        #undo_relax=JuMP.relax_integrality(pm.model)
        if (haskey(pm.setting,"agent") && pm.setting["agent"]!="")
            fix_dc_ne_lines2zero(pm)
            fix_ac_ne_lines2zero(pm)
            fix_variables(pm)
            objective_min_cost_acdc_convexcble_strg_npv_admm(pm)
        else
            if (haskey(pm.setting,"relax_problem") && pm.setting["relax_problem"]==true)
                objective_min_cost_acdc_convex_convcble_strg_npv(pm)
                fix_dc_ne_lines2zero(pm)
                fix_ac_ne_lines2zero(pm)
            else
                objective_min_cost_acdc_convex_conv_strg_npv(pm)
            end
        end
        #OBJECTIVE see objective.jl

        #OBJECTIVE see objective.jl
        #objective_min_cost_acdc_convex_conv_strg_npv(pm)
        #CONSTRAINTS: defined within PowerModels(ACDC) can directly be used, other constraints need to be defined in the according sections of the code: flexible_demand.jl
        #=for v in JuMP.all_variables(pm.model)
            println(v)
        end=#
        #println(JuMP.objective_function(pm.model))
        #println(pm.model)

end

function admm_4_AjAwAgAuAo_main(mn_data, gurobi, s)
    results_set=[]
    eps=s["eps"]#set max value of residual for convergence
    residual=Inf
    agents=["Ag","Au","Ao","Aw","Aj"]#Initilize agants
    push!(s,"fixed_variables" => Dict{String,Any}())
    push!(s,"agent" => "")
    s["fixed_variables"] = admm_4_AjAwAgAuAo_intialize(mn_data["nw"], s["fixed_variables"], s["wfz"], s["genz"])
    while (residual>eps)
    #for i=1:10
        results_set=[]
        for a in agents
            s["agent"]=a
            println(a)
            result_mip = cordoba_acdc_wf_strg_admm(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
            if !(isnan(result_mip["objective"]))
                s["fixed_variables"] = admm_4_AjAwAgAuAo_update(result_mip["solution"]["nw"], s["fixed_variables"], a, s["wfz"])
                push!(results_set,(a,result_mip))
            else
                println("skyap!!!")
            end
        end
        [println(first(r)*" "*string(last(r)["objective"])) for r in results_set]
        s["fixed_variables"], residual = dual_variable_update(s["fixed_variables"])
        println("Residual: "*string(residual))
    #end
    end
    return last(results_set[length(results_set)])
end


function fix_dc_lines2zero(pm)
    for n in _PM.nw_ids(pm)
        if haskey(_PM.ref(pm, n), :branchdc)
            branchdc = _PM.ref(pm, n, :branchdc)
            if !isempty(branchdc)
                for (i,br) in branchdc
                    brnch=_PM.var(pm, n, :branchdc, i)
                    JuMP.fix(brnch,0,force=true)
                end
            end
        end
    end
end

function fix_ac_lines2zero(pm)
    for n in _PM.nw_ids(pm)
        if haskey(_PM.ref(pm, n), :branch)
            branch = _PM.ref(pm, n, :branch)
            if !isempty(branch)
                for (i,br) in branch
                    brnch=_PM.var(pm, n, :branch, i)
                    JuMP.fix(brnch,0,force=true)
                end
            end
        end
    end
end


function fix_dc_ne_lines2zero(pm)
    for n in _PM.nw_ids(pm)
        if haskey(_PM.ref(pm, n), :branchdc_ne)
            branchdc_ne = _PM.ref(pm, n, :branchdc_ne)
            if !isempty(branchdc_ne)
                for (i,br) in branchdc_ne
                    brnch=_PM.var(pm, n, :branchdc_ne, i)
                    JuMP.fix(brnch,0,force=true)
                end
            end
        end
    end
end

function fix_ac_ne_lines2zero(pm)
    for n in _PM.nw_ids(pm)
        if haskey(_PM.ref(pm, n), :branch_ne)
            branch_ne = _PM.ref(pm, n, :branch_ne)
            if !isempty(branch_ne)
                for (i,br) in branch_ne
                    brnch=_PM.var(pm, n, :branch_ne, i)
                    JuMP.fix(brnch,0,force=true)
                end
            end
        end
    end
end

#=function constraint_ohms_yt_from(pm::_PM.AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = _PM.var(pm, n,  :p, f_idx)
    va_fr = _PM.var(pm, n, :va, f_bus)
    va_to = _PM.var(pm, n, :va, t_bus)

    JuMP.@constraint(pm.model, p_fr == -b*(va_fr - va_to))
    # omit reactive constraint
end

function constraint_ohms_yt_to(pm::_PM.AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_to  = _PM.var(pm, n,  :p, t_idx)
    va_fr = _PM.var(pm, n, :va, f_bus)
    va_to = _PM.var(pm, n, :va, t_bus)

    JuMP.@constraint(pm.model, p_to == -b*(va_to - va_fr))
    # omit reactive constraint
end=#
function constraint_ohms_dc_branch_ne(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = PowerModels.ref(pm, nw, :branchdc_ne, i)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    rate_a = branch["rateA"]
    #println(rate_a)
    #println("rateA")
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    p = PowerModels.ref(pm, nw, :dcpol)
#    if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
#        println("%%%%% Inside!!!!!!!!!! %%%%%")
    #    println(pm.setting["home_market"])
        constraint_ohms_dc_branch_ne(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["r"], p, rate_a)
#    end
end

function constraint_ohms_dc_branch_ne(pm::_PM.AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, r, p, rate_a)
    p_dc_fr_ne = _PM.var(pm, n, :p_dcgrid_ne, f_idx)
    p_dc_to_ne = _PM.var(pm, n, :p_dcgrid_ne, t_idx)
#    if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
        JuMP.@constraint(pm.model, p_dc_fr_ne + p_dc_to_ne == 0)
#    else
#        JuMP.@constraint(pm.model, p_dc_fr_ne + p_dc_to_ne + 0.3*rate_a >= 0)
#        JuMP.@constraint(pm.model, p_dc_fr_ne + p_dc_to_ne - 0.3*rate_a <= 0)
#    end
end

function constraint_ohms_dc_branch(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PM.ref(pm, nw, :branchdc, i)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    p_rateA=_PM.var(pm, nw, :p_rateA, i)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)
    #if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
        constraint_ohms_dc_branch(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["r"], p, p_rateA)
    #end
end

function constraint_ohms_dc_branch(pm::_PM.AbstractDCPModel, n::Int,  f_bus, t_bus, f_idx, t_idx, r, p, p_rateA)
    p_dc_fr = _PM.var(pm, n, :p_dcgrid, f_idx)
    p_dc_to = _PM.var(pm, n, :p_dcgrid, t_idx)

    #if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
        JuMP.@constraint(pm.model, p_dc_fr + p_dc_to == 0)
    #else
    #    println("p_rateA")
    #    println(p_rateA)
    #    JuMP.@constraint(pm.model, p_dc_fr + p_dc_to + 0.3*p_rateA >= 0)
    #    JuMP.@constraint(pm.model, p_dc_fr + p_dc_to - 0.3*p_rateA <= 0)
    #end
end

function constraint_ohms_ac_branch(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PM.ref(pm, nw, :branch, i)
    rate_a =_PM.var(pm, nw, :p_rateAC, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_ohms_ac_branch(pm, nw, f_idx, t_idx, f_bus, t_bus, rate_a)
end

function constraint_ohms_ac_branch(pm::_PM.AbstractDCPModel, n::Int, f_idx, t_idx, f_bus, t_bus, rate_a)
    p_fr  = _PM.var(pm, n,  :p, f_idx)
    p_to  = _PM.var(pm, n,  :p, t_idx)
#    if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
        JuMP.@constraint(pm.model, p_fr + p_to == 0)
#    else
#        JuMP.@constraint(pm.model, p_fr + p_to + 0.3*rate_a >= 0)
#        JuMP.@constraint(pm.model, p_fr + p_to - 0.3*rate_a <= 0)
#    end
end

function constraint_ohms_ac_branch_ne(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    rate_a = branch["rate_a"]
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    #if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
        constraint_ohms_ac_branch_ne(pm, nw, f_idx, t_idx, f_bus, t_bus, rate_a)
    #end
end

function constraint_ohms_ac_branch_ne(pm::_PM.AbstractDCPModel, n::Int, f_idx, t_idx, f_bus, t_bus, rate_a)
    #z = _PM.var(pm, n, :branch_ne, i)
    p_fr  = _PM.var(pm, n,  :p_ne, f_idx)
    p_to  = _PM.var(pm, n,  :p_ne, t_idx)
#    if (!(haskey(pm.setting, "home_market")) || !(issubset([f_bus],pm.setting["home_market"]) && issubset([t_bus],pm.setting["home_market"])))
        JuMP.@constraint(pm.model, p_fr + p_to == 0)
#    else
#        JuMP.@constraint(pm.model, p_fr + p_to + 0.3*rate_a >= 0)
#        JuMP.@constraint(pm.model, p_fr + p_to - 0.3*rate_a <= 0)
#    end
end

#=
function constraint_ohms_yt_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    _PM.constraint_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_ohms_yt_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    _PM.constraint_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end
=#
