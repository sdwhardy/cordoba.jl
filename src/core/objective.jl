##################################################################
##################### Objective with candidate storage
##################################################################

##################################### NPV considered #############################
#binary AC and DC transmission lines, continuous converter (NPV considered)
function objective_min_cost_acdc_convex_conv_npv(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(pm.ref[:scenario_prob][s] *
            sum(
                calc_gen_cost(pm, n)
                + calc_convdc_convexafy_cost_npv(pm, n)
                + calc_ne_branch_cost(pm, n)#was commented out - not tested since uncommneted
                + calc_branchdc_ne_cost(pm, n)
            for (sc, n) in scenario)
        for (s, scenario) in pm.ref[:scenario])
    )
end

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

#Constraint on size of yearly investment possible
#=function max_investment_per_year_old(pm::_PM.AbstractPowerModel)
    s1=pm.setting["scenarios_length"]
    y1=pm.setting["years_length"]
    h1=pm.setting["hours_length"]
    #println(pm.ref[:scenario])
    for (s, scenario) in pm.ref[:scenario]
        s_num=parse(Int64,s)
        jan=(s_num-1)*y1*h1+1
        while (jan<=(s_num)*y1*h1)
            cost=0.0
            dec=jan+h1-1
            for (sc, n) in scenario
                if (jan<=n && n<=dec)
                    cost=cost+calc_convdc_convexafy_cost_npv(pm, n)
                    cost=cost+calc_ne_branch_cost(pm, n)
                    cost=cost+calc_branchdc_ne_cost(pm, n)
                    cost=cost+calc_storage_cost_cordoba_npv(pm, n)
                    cost=cost+calc_wf_cost_npv(pm, n)
                end
            end
            JuMP.@constraint(pm.model,cost<=pm.setting["max_invest_per_year"])
            jan=jan+h1
        end;
    end
end=#

#Constraint on size of yearly investment possible
function max_investment_per_year(pm::_PM.AbstractPowerModel)
    s1=pm.setting["scenarios_length"]
    y1=pm.setting["years_length"]
    h1=pm.setting["hours_length"]
    #println(pm.ref[:scenario])
    for (s, scenario) in pm.ref[:scenario]
        s_num=parse(Int64,s)
        #jan=(s_num-1)*y1*h1+1
        #while (jan<=(s_num)*y1*h1)
            cost=0.0
            #dec=jan+h1-1
            ordered_scenes=sort(OrderedCollections.OrderedDict(scenario), by=x->parse(Int64,x))
            ks=collect(keys(ordered_scenes))
            strt=ordered_scenes[ks[1]]
            _sc=floor(Int64,(strt-1)/(y1*h1))
            _yr=ceil(Int64,(strt-_sc*(y1*h1))/(h1))
            for (sc, n) in ordered_scenes
                if (pm.setting["relax_problem"])
                    cost=cost+calc_branchdc_cost_max_invest(pm, n)
                    if (pm.setting["AC"]=="1")
                        cost=cost+calc_branch_cost_max_invest(pm, n);end
                else
                    cost=cost+calc_branchdc_ne_cost_max_invest(pm, n)
                    if (pm.setting["AC"]=="1")
                        cost=cost+calc_branch_ne_cost_max_invest(pm, n);
                    end
                end
                cost=cost+calc_convdc_convexafy_cost_max_invest(pm, n)
                cost=cost+calc_storage_cost_cordoba_max_invest(pm, n)
                cost=cost+calc_wf_cost_max_invest(pm, n)
                if (n==strt+h1-1)
                    JuMP.@constraint(pm.model,cost<=pm.setting["max_invest_per_year"][_yr])
                    strt=strt+h1
                    _yr=ceil(Int64,(strt-_sc*(y1*h1))/(h1))
                    cost=0.0
                end
                #println("n: "*string(n))
            end

            #jan=jan+h1
        #end;
    end
end

function calc_branch_ne_cost_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branch_ne_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:branch_ne,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    brs=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(brs,_PM.ref(pm, nt, :ne_branch));end
        for bs in brs
            cost = cost + sum(calc_single_branch_ne_cost_npv(i,b["construction_cost"],n) for (i,b) in bs)
        end
    else
        for nt=n:hl:n-hl+(yl-_yr+1)*hl
            push!(brs,_PM.ref(pm, nt, :ne_branch));end
        for bs in brs
            cost = cost + sum(calc_single_branch_ne_cost_npv(i,b["construction_cost"],n) for (i,b) in bs)
            cost = cost - sum(calc_single_branch_ne_cost_npv(i,b["construction_cost"],n-hl) for (i,b) in bs)
        end
    end
    return cost
end

function calc_branchdc_ne_cost_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branchdc_ne_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:branchdc_ne,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    brs=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(brs,_PM.ref(pm, nt, :branchdc_ne));end
        for bs in brs
            cost = cost + sum(calc_single_branchdc_ne_cost_npv(i,b["cost"],n) for (i,b) in bs)
        end
    else
        for nt=n:hl:n-hl+(yl-_yr+1)*hl
            push!(brs,_PM.ref(pm, nt, :branchdc_ne));end
        for bs in brs
            cost = cost + sum(calc_single_branchdc_ne_cost_npv(i,b["cost"],n) for (i,b) in bs)
            cost = cost - sum(calc_single_branchdc_ne_cost_npv(i,b["cost"],n-hl) for (i,b) in bs)
        end
    end
    return cost
end

function calc_branchdc_cost_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branchdc_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:p_rateA,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    brs=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(brs,_PM.ref(pm, nt, :branchdc));end
        for bs in brs
            cost = cost + sum(calc_single_branchdc_cost_npv(i,b["cost"],n) for (i,b) in bs)
        end
    else
        for nt=n:hl:n+(yl-_yr)*hl
            push!(brs,_PM.ref(pm, nt, :branchdc));end
        for bs in brs
            cost = cost + sum(calc_single_branchdc_cost_npv(i,b["cost"],n) for (i,b) in bs)
            cost = cost - sum(calc_single_branchdc_cost_npv(i,b["cost"],n-hl) for (i,b) in bs)
        end
    end
    return cost
end

function calc_branch_cost_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branch_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:p_rateAC,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    brs=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(brs,_PM.ref(pm, nt, :branch));end
        for bs in brs
            cost = cost + sum(calc_single_branch_cost_npv(i,b["cost"],n) for (i,b) in bs)
        end
    else
        for nt=n:hl:n+(yl-_yr)*hl
            push!(brs,_PM.ref(pm, nt, :branch));end
        for bs in brs
            cost = cost + sum(calc_single_branch_cost_npv(i,b["cost"],n) for (i,b) in bs)
            cost = cost - sum(calc_single_branch_cost_npv(i,b["cost"],n-hl) for (i,b) in bs)
        end
    end
    return cost
end

#max investment per year for converters
function calc_wf_cost_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_wf_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:wf_pacmax,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    wfs=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(wfs,_PM.ref(pm, nt, :gen));end
        for gs in wfs
            cost = cost + sum(calc_single_wf_cost_npv(i,g["invest"],n) for (i,g) in gs if issubset([i],first.(pm.setting["wfz"])))
        end
    else
        for nt=n:hl:n+(yl-_yr)*hl
            push!(wfs,_PM.ref(pm, nt, :gen));end
        for gs in wfs
            cost = cost + sum(calc_single_wf_cost_npv(i,g["invest"],n) for (i,g) in gs if issubset([i],first.(pm.setting["wfz"])))
            cost = cost - sum(calc_single_wf_cost_npv(i,g["invest"],n-hl) for (i,g) in gs if issubset([i],first.(pm.setting["wfz"])))
        end
    end
    return cost
end

#max investment per year for converters
function calc_storage_cost_cordoba_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_storage_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:e_absmax,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    stores=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(stores,_PM.ref(pm, nt, :storage));end
        for s in stores
            cost = cost + sum(calc_single_storage_cost_npv(i,b["cost"],n) for (i,b) in s)
        end
    else
        for nt=n:hl:n+(yl-_yr)*hl
            push!(stores,_PM.ref(pm, nt, :storage));end
        for s in stores
            cost = cost + sum(calc_single_storage_cost_npv(i,b["cost"],n) for (i,b) in s)
            cost = cost - sum(calc_single_storage_cost_npv(i,b["cost"],n-hl) for (i,b) in s)
        end
    end
    return cost
end

#max investment per year for converters
function calc_convdc_convexafy_cost_max_invest(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_convdc_cost_npv(i, b_cost, nw)
        cost = 0.0
        #println(string(nw)*" "*string(i)*" "*string(b_cost))
        cost += b_cost * _PM.var(pm,nw,:p_pacmax,i)
        return cost
    end

    sl=pm.setting["scenarios_length"]
    yl=pm.setting["years_length"]
    hl=pm.setting["hours_length"]
    _sc=floor(Int64,(n-1)/(yl*hl))
    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
    convdc=[]
    cost = 0
    if (_yr==1)
        for nt=n:hl:n+yl*hl-1
            push!(convdc,_PM.ref(pm, nt, :convdc));end
        for c in convdc
            cost = cost + sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in c)
        end
    else
        for nt=n:hl:n+(yl-_yr)*hl
            push!(convdc,_PM.ref(pm, nt, :convdc));end
        for c in convdc
            cost = cost + sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in c)
            cost = cost - sum(calc_single_convdc_cost_npv(i,b["cost"],n-hl) for (i,b) in c)
        end
    end
    return cost
end

#convex converter considering NPV
function calc_convdc_convexafy_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_convdc_cost_npv(i, b_cost, nw)
        cost = 0.0
        cost += b_cost * _PM.var(pm,nw,:p_pacmax,i)
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
        cost = cost-sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in convdc1)
    else
        cost = sum(calc_single_convdc_cost_npv(i,b["cost"],n) for (i,b) in convdc0)
    end
    return cost
end


#cost of expanding storage year by year NPV
function calc_storage_cost_cordoba_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_storage_cost_npv(i, b_cost, nw)
        cost = 0.0
        cost += b_cost * _PM.var(pm,nw,:e_absmax,i)
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
        cost = 0.0
        cost += b_cost * _PM.var(pm,nw,:wf_pacmax,i)
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


################################## Common ######################################################
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

################################ No NPV considered ######################################
#Only converters and dc lines considered - not considering NPV
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

#convex converter not considering NPV
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
#convex ac branch
function calc_branch_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_branch_cost_npv(i, b_cost, nw)
        cost = 0.0
        cost += b_cost * _PM.var(pm,nw,:p_rateAC,i) #
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
        cost = 0.0
        cost += b_cost * _PM.var(pm,nw,:p_rateA,i) #
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

#convex dc branch - No NPV
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

#convex converter considering NPV
function calc_convdc_convexafy_cost_npv(pm::_PM.AbstractPowerModel, n::Int)

    function calc_single_convdc_cost_npv(i, b_cost, nw)
        cost = 0.0
        cost += b_cost * _PM.var(pm,nw,:p_pacmax,i)
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
##################################### Binary candidates ########################
#all candidates are binaries
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

#binary storage with inj/abs costs included - for FFR analysis
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

####################################### Depricated #############################

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

#=function objective_min_cost_acdc_cordoba(pm::_PM.AbstractPowerModel)
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
end=#


#=function objective_min_cost_acdc_convex_conv(pm::_PM.AbstractPowerModel)
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
end=#

#=function calc_storage_cost_cordoba_npv(pm::_PM.AbstractPowerModel, n::Int)
    #addition of buy and sell price storage["cost_abs"],storage["cost_inj"]
    function calc_single_storage_cost(i, cost_a,cost_i)
        #if (_PM.var(pm,n,:ps,i)>0)
            cost = -1*cost_i * _PM.var(pm, n, :sd_ne, i)*_PM.var(pm, n, :z_strg_ne, i)# takes profit
        #else
            cost += -1*cost_a * _PM.var(pm, n, :sc_ne, i)*_PM.var(pm, n, :z_strg_ne, i)# buys energy
        #end
        return cost
    end

    if haskey(_PM.ref(pm, n), :storage)
        storage = _PM.ref(pm, n, :storage)

        if !isempty(storage)
            #cost for buy and sell MWh
            cost = sum(calc_single_storage_cost(i,strg["cost"]) for (i,strg) in storage)
            cost += sum((storage["eq_cost"] + storage["inst_cost"])*_PM.var(pm, n, :z_strg_ne, i) for (i,storage) in ne_storage)
            if add_co2_cost
                cost += sum(storage["co2_cost"]*_PM.var(pm, n, :z_strg_ne, i) for (i,storage) in ne_storage)
            end
        end
    end
    return cost
end=#
