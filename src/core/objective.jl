##################################################################
##################### Objective with candidate storage
##################################################################

function objective_min_cost_storage_cordoba(pm::_PM.AbstractPowerModel)
    add_co2_cost = haskey(pm.setting, "add_co2_cost") && pm.setting["add_co2_cost"]
    return JuMP.@objective(pm.model, Min,
        sum(
            _FP.calc_gen_cost(pm, n, add_co2_cost)
            + _FP.calc_convdc_ne_cost(pm, n, add_co2_cost)
            + _FP.calc_ne_branch_cost(pm, n, add_co2_cost)
            + _FP.calc_branchdc_ne_cost(pm, n, add_co2_cost)
            + calc_ne_storage_cost_cordoba(pm, n, add_co2_cost)
        for n in _PM.nw_ids(pm))
    )
end

function calc_ne_storage_cost_cordoba(pm::_PM.AbstractPowerModel, n::Int, add_co2_cost::Bool)
    #addition of buy and sell price storage["cost_abs"],storage["cost_inj"]
    function calc_single_ne_storage_cost(i, cost_a,cost_i)
        #if (_PM.var(pm,n,:ps,i)>0)
            cost = -1*cost_i * _PM.var(pm, n, :sd_ne, i)*_PM.var(pm, n, :z_strg_ne, i)# takes profit
        #else
            cost += -1*cost_a * _PM.var(pm, n, :sc_ne, i)*_PM.var(pm, n, :z_strg_ne, i)# buys energy
        #end
        return cost
    end

    if haskey(_PM.ref(pm, n), :ne_storage)
        ne_storage = _PM.ref(pm, n, :ne_storage)

        if !isempty(ne_storage)
            #cost for buy and sell MWh
            cost = sum(calc_single_ne_storage_cost(i,storage["cost_abs"],storage["cost_inj"]) for (i,storage) in ne_storage)
            cost += sum((storage["eq_cost"] + storage["inst_cost"])*_PM.var(pm, n, :z_strg_ne, i) for (i,storage) in ne_storage)
            if add_co2_cost
                cost += sum(storage["co2_cost"]*_PM.var(pm, n, :z_strg_ne, i) for (i,storage) in ne_storage)
            end
        end
    end
    return cost
end
