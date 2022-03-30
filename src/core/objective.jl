#Complete cordoba objective function
#binary AC and DC transmission lines, continuous converter, continuous storage, wind farm expansion (NPV considered)
function objective_min_cost_acdc_convex_conv_strg_npv(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                calc_gen_cost(pm, n)
                + calc_convdc_convexafy_cost_npv(pm, n)
                + calc_ne_branch_cost(pm, n)
                + calc_branchdc_ne_cost(pm, n)
                + calc_storage_cost_cordoba_npv(pm, n)
                + calc_wf_cost_npv(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

#Objective with convex Cables, continuous converter, continuous storage, wind farm expansion (NPV considered) 
function objective_min_cost_acdc_convex_convcble_strg_npv(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                calc_gen_cost(pm, n)
                + calc_convdc_convexafy_cost_npv(pm, n)
                + calc_branch_cost_npv(pm, n)
                + calc_branchdc_cost_npv(pm, n)
                + calc_storage_cost_cordoba_npv(pm, n)
                + calc_wf_cost_npv(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

################################### Convex variable costs #######################################
#cost of generation
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

#convex converter considering NPV
function calc_convdc_convexafy_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_convdc_cost_npv(i, b_cost, nw)
        cost = b_cost * _PM.var(pm,nw,:p_pacmax,i)
        return cost
    end
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    convdc0 = _PM.ref(pm, n, :convdc)

    if (_yr<yl)
        #println("if "*string(n)*" "*string(hl)*" "*string(yl)*" "*string(_yr))
        convdc1 = _PM.ref(pm, n+hl, :convdc)
        cost = sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in convdc0)
        #println(cost)
        cost = cost-sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in convdc1)
        #println(cost)
    else
        #println("else "*string(n)*" "*string(hl)*" "*string(yl)*" "*string(_yr))
        cost = sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in convdc0)
    end
    return cost
end


#cost of expanding storage year by year NPV
function calc_storage_cost_cordoba_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_storage_cost_npv(i, b_cost, nw)
        cost = b_cost * _PM.var(pm,nw,:e_absmax,i)       
        return cost
    end
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    strg0 = _PM.ref(pm, n, :storage)

    if (_yr<yl)
        strg1 = _PM.ref(pm, n+hl, :storage)
        cost = sum(calc_single_storage_cost_npv(i,s["cost"],n) for (i,s) in strg0)
        cost = cost+sum(calc_single_storage_cost_npv(i,s["cost"]*-1,n) for (i,s) in strg1)
    else
        cost = sum(calc_single_storage_cost_npv(i,s["cost"],n) for (i,s) in strg0)
    end
    return cost
end

#cost of expanding a wind farm year by year NPV
function calc_wf_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_wf_cost_npv(i, b_cost,nw)       
        cost = b_cost * _PM.var(pm,nw,:wf_pacmax,i)
        return cost
    end
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    gen0 = _PM.ref(pm, n, :gen)

    if (_yr<yl)
        gen1 = _PM.ref(pm, n+hl, :gen)
        cost = sum(calc_single_wf_cost_npv(i,g["invest"],n) for (i,g) in gen0 if issubset([i],first.(pm.setting["wfz"])))
        cost = cost+sum(calc_single_wf_cost_npv(i,g["invest"]*-1,n) for (i,g) in gen1 if issubset([i],first.(pm.setting["wfz"])))
    else
        cost = sum(calc_single_wf_cost_npv(i,g["invest"],n) for (i,g) in gen0 if issubset([i],first.(pm.setting["wfz"])))
    end
    return cost
end

#convex ac branch
function calc_branch_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branch_cost_npv(i, b_cost, nw)  
        cost = b_cost * _PM.var(pm,nw,:p_rateAC,i) #
        return cost
    end
    if (pm.setting["AC"]=="1")
        sl=pm.setting["scenarios_length"]
        yl=pm.setting["years_length"]
        hl=pm.setting["hours_length"]
        _sc=floor(Int64,(n-1)/(yl*hl))
        _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
        branchdc0 = _PM.ref(pm, n, :branch)
        if (_yr<yl)
            branchdc1 = _PM.ref(pm, n+hl, :branch)
            cost = sum(calc_single_branch_cost_npv(i,b["cost"],n) for (i,b) in branchdc0)
            cost = cost-sum(calc_single_branch_cost_npv(i,b["cost"],n) for (i,b) in branchdc1)
        else
            cost = sum(calc_single_branch_cost_npv(i,b["cost"],n) for (i,b) in branchdc0)
        end
    else
        cost=0
    end
    return cost
end

#convex dc branch
function calc_branchdc_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branchdc_cost_npv(i, b_cost, nw)
        cost = b_cost * _PM.var(pm,nw,:p_rateA,i) 
        return cost
    end
    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    branchdc0 = _PM.ref(pm, n, :branchdc)
    if (_yr<yl)
        branchdc1 = _PM.ref(pm, n+hl, :branchdc)
        cost = sum(calc_single_branchdc_cost_npv(i,b["cost"],n) for (i,b) in branchdc0)
        cost = cost-sum(calc_single_branchdc_cost_npv(i,b["cost"],n) for (i,b) in branchdc1)
    else
        cost = sum(calc_single_branchdc_cost_npv(i,b["cost"],n) for (i,b) in branchdc0)
    end
    return cost
end

################################## Binary variable costs ###############################
#cost of an AC branch
function calc_ne_branch_cost(pm::_PM.AbstractPowerModel, n::Int)
    cost = 0.0
    if haskey(_PM.ref(pm, n), :ne_branch)
        ne_branch = _PM.ref(pm, n, :ne_branch)
        if !isempty(ne_branch)
            #println(_PM.var(pm, n, :branch_ne, i) for (i,branch) in ne_branch)
            cost = sum(branch["construction_cost"]*_PM.var(pm, n, :branch_ne, i) for (i,branch) in ne_branch)
            #println(cost)
        end
    end
    return cost
end

#cost of a dc branch
function calc_branchdc_ne_cost(pm::_PM.AbstractPowerModel, n::Int)
    cost = 0.0
    if haskey(_PM.ref(pm, n), :branchdc_ne)
        branchdc_ne = _PM.ref(pm, n, :branchdc_ne)
        if !isempty(branchdc_ne)
            #println(_PM.var(pm, n, :branchdc_ne, i) for (i,branch) in branchdc_ne)
            cost = sum(branch["cost"]*_PM.var(pm, n, :branchdc_ne, i) for (i,branch) in branchdc_ne)
            #println(cost)
        end
    end
    return cost
end
