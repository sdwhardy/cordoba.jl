##################################################################
##################### Objective with candidate storage
##################################################################
function objective_min_cost_acdc_cordoba(pm::_PM.AbstractPowerModel)
    add_co2_cost = haskey(pm.setting, "add_co2_cost") && pm.setting["add_co2_cost"]
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                _FP.calc_gen_cost(pm, n, add_co2_cost)
                + _FP.calc_convdc_ne_cost(pm, n, add_co2_cost)
                + _FP.calc_ne_branch_cost(pm, n, add_co2_cost)
                + _FP.calc_branchdc_ne_cost(pm, n, add_co2_cost)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

function objective_min_cost_opf(pm::_PM.AbstractPowerModel)
    add_co2_cost = haskey(pm.setting, "add_co2_cost") && pm.setting["add_co2_cost"]
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                _FP.calc_gen_cost(pm, n, add_co2_cost)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

function objective_min_cost_acdc_convex_conv(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                calc_gen_cost(pm, n)
                + calc_convdc_convexafy_cost(pm, n)
                + calc_ne_branch_cost(pm, n)
                + calc_branchdc_ne_cost(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

function objective_min_cost_acdc_convex_conv_npv(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                #calc_gen_cost(pm, n)
                calc_convdc_convexafy_cost_npv(pm, n)
                #+ calc_ne_branch_cost(pm, n)
                #+ calc_branchdc_ne_cost(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

function calc_branchdc_ne_cost(pm::_PM.AbstractPowerModel, n::Int)
    cost = 0.0
    if haskey(_PM.ref(pm, n), :branchdc_ne)
        branchdc_ne = _PM.ref(pm, n, :branchdc_ne)
        if !isempty(branchdc_ne)
            cost = sum(branch["cost"]*_PM.var(pm, n, :branchdc_ne, i) for (i,branch) in branchdc_ne)
            #println(cost)
        end
    end
    return cost
end

function calc_ne_branch_cost(pm::_PM.AbstractPowerModel, n::Int)
    cost = 0.0
    if haskey(_PM.ref(pm, n), :ne_branch)
        ne_branch = _PM.ref(pm, n, :ne_branch)
        if !isempty(ne_branch)
            cost = sum(branch["construction_cost"]*_PM.var(pm, n, :branch_ne, i) for (i,branch) in ne_branch)
        end
    end
    return cost
end

function calc_gen_cost(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_gen_cost(i, g_cost)
        len = length(g_cost)
        cost = 0.0
        if len >= 1
            cost = g_cost[len] # Constant term
            if len >= 2
                cost += g_cost[len-1] * _PM.var(pm,n,:pg,i) # Adds linear term
            end
        end
        return cost
    end

    gen = _PM.ref(pm, n, :gen)
    cost = sum(calc_single_gen_cost(i,g["cost"]) for (i,g) in gen)
    return cost
end


function objective_min_cost_acdc_convexafy(pm::_PM.AbstractPowerModel)
    add_co2_cost = haskey(pm.setting, "add_co2_cost") && pm.setting["add_co2_cost"]
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                _FP.calc_gen_cost(pm, n, add_co2_cost)
                + calc_convdc_convexafy_cost(pm, n)
                + calc_branchdc_ne_cost(pm, n)
                #+ calc_branchdc_convexafy_cost(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

function calc_convdc_convexafy_cost(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_convdc_cost(i, b_cost)
        cost = 0.0
        cost += b_cost * _PM.var(pm,n,:p_pacmax,i) #
        return cost
    end

    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    convdc = _PM.ref(pm, n, :convdc)
    if (mod(n-1,yl*hl)>=(yl-1)*hl)
        #println("n: "*string(n)*", mod(n,yl*hl): "*string(mod(n,yl*hl))*", hl: "*string(hl))
        cost = sum(calc_single_convdc_cost(i,b["cost"]) for (i,b) in convdc)
    else
        cost = 0#sum(calc_single_branchdc_yr2on_cost(i,b["rateC"],hl) for (i,b) in branchdc)
    end
    return cost
end




function calc_convdc_convexafy_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_convdc_cost_npv(i, b_cost)
        cost = 0.0
        cost += b_cost * _PM.var(pm,n,:p_pacmax,i)
        return cost
    end
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    convdc0 = _PM.ref(pm, n, :convdc)
    #=k0=keys(convdc0)
    c0=[k for k in k0]
    println(string(n)*" convdc0: "*string(convdc0[c0[1]]["cost"])*" - "*string(convdc0[c0[2]]["cost"]))=#
    if (_yr<yl)
        convdc1 = _PM.ref(pm, n+hl, :convdc)
        #=k1=keys(convdc1)
        c1=[k for k in k1]
        println(string(n+hl)*" convdc1: "*string(convdc1[c1[1]]["cost"])*" - "*string(convdc1[c1[2]]["cost"]))=#
        cost = sum(calc_single_convdc_cost_npv(i,b["cost"]) for (i,b) in convdc0)
        cost = cost+sum(calc_single_convdc_cost_npv(i,b["cost"]*-1) for (i,b) in convdc1)
    else
        cost = sum(calc_single_convdc_cost_npv(i,b["cost"]) for (i,b) in convdc0)
    end
    return cost
end

function calc_branchdc_convexafy_cost(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branchdc_cost(i, b_cost)
        cost = 0.0
        cost += b_cost * _PM.var(pm,n,:p_rateA,i) #
        return cost
    end
    function calc_single_branchdc_yr2on_cost(i, b_cost,hl)
        cost = 0.0
        cost += b_cost * _PM.var(pm,n,:p_rateA,i)
        cost -= b_cost * _PM.var(pm,n-hl,:p_rateA,i) #
        return cost
    end
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    branchdc = _PM.ref(pm, n, :branchdc)
    if (mod(n-1,yl*hl)>=(yl-1)*hl)
        #println("n: "*string(n)*", mod(n,yl*hl): "*string(mod(n,yl*hl))*", hl: "*string(hl))
        cost = sum(calc_single_branchdc_cost(i,b["cost"]) for (i,b) in branchdc)
    else
        cost = 0#sum(calc_single_branchdc_yr2on_cost(i,b["rateC"],hl) for (i,b) in branchdc)
    end
    return cost
end

function calc_branchdc_cost(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branchdc_cost(i, b_cost)
        cost = 0.0
        cost += b_cost * _PM.var(pm,n,:p_rateA,i) #
        return cost
    end

    branchdc = _PM.ref(pm, n, :branchdc)
    cost = sum(calc_single_branchdc_cost(i,b["cost"]) for (i,b) in branchdc)
    return cost
end

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
