
function objective_min_cost_acdc_convexcble_strg_npv_admm(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                calc_gen_cost(pm, n)
                + calc_admm_penalty_cost(pm, n)
                + calc_convdc_convexafy_cost_npv(pm, n)
                + calc_branch_cost_npv(pm, n)
                + calc_branchdc_cost_npv(pm, n)
                + calc_storage_cost_cordoba_npv(pm, n)
                + calc_wf_cost_npv(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

function objective_min_cost_acdc_convex_conv_strg_npv_admm(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                calc_gen_cost(pm, n)
                + calc_admm_penalty_cost(pm, n)
                + calc_convdc_convexafy_cost_npv(pm, n)
                + calc_ne_branch_cost(pm, n)
                + calc_branchdc_ne_cost(pm, n)
                + calc_storage_cost_cordoba_npv(pm, n)
                + calc_wf_cost_npv(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

#cost of generation
function calc_admm_penalty_cost(pm::_PM.AbstractPowerModel, n::Int)
    cost = 0.0
    for i in _PM.ids(pm, n, :bus)
        if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
            if !(issubset(i,pm.setting["home_market"]))
                cost += constraint_power_balance_acne_dcne_strg(pm, i; nw = n)
            end
        else
            cost += constraint_power_balance_acne_dcne_strg(pm, i; nw = n)
        end
    end
    if haskey(pm.setting, "home_market") && length(pm.setting["home_market"]) > 1
        is=intersect(_PM.ids(pm, n, :bus),first.(pm.setting["home_market"]))
        cost += constraint_power_balance_acne_dcne_strg_hm(pm, is; nw = n)
    end
    return cost
end
