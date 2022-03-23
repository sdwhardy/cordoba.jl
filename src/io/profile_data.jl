#load Time series data
function load_time_series(rt_ex, argz)
    scenario_data=FileIO.load(rt_ex*"time_series_k"*string(argz["k"])*".jld2")
    #keep only specified scenarios
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],argz["scenario_names"]));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],argz["scenario_years"]));delete!(scenario_data[k],y);end; end;end;end
    if (haskey(argz, "test") && argz["test"]==true)
        for k0 in d_keys; for k1 in keys(scenario_data[k0]); scenario_data[k0][k1]=scenario_data[k0][k1][1:2,:];end;end
    end
    ##################### Find minimum length scenario and Make all scenarios the same length
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end
    return scenario_data, ls
end


#multi period problem setup
function multi_period_setup(ls,scenario_data,data, markets, infinite_grid, argz)
    #################### Multi-period input parameters #######################
    all_scenario_data,data,scenario, dim = multi_period_stoch_year_setup(ls,argz["scenario_years"],argz["scenario_names"],scenario_data,data)
    scenario["planning_horizon"] = argz["scenario_planning_horizon"] # in years, to scale generation cost
    extradata,data =create_profile_sets_mesh(dim, data, all_scenario_data, markets, infinite_grid, [data["baseMVA"] for wf in argz["owpp_mva"]])
    #########################################################################
    #################### Scale cost data
    #[println(b*" "*string(br["construction_cost"])) for (b,br) in data["ne_branch"]];println()
    scale_cost_data_2hourly!(data, scenario)#infrastructure investments
    #[println(b*" "*string(br["construction_cost"])) for (b,br) in data["ne_branch"]];println()
    scale_cost_data_2yearly!(extradata, scenario)#energy cost benefits
    #[println(b*" "*string(br["construction_cost"])) for (b,br) in data["ne_branch"]];println()

    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    #[println(k*" "*b*" "*string(br["construction_cost"])) for (k,nw) in mn_data["nw"] for (b,br) in nw["ne_branch"]];println()
    # scale all to NPV
    #mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)#sum of years must equal total
    mn_data = npvs_costs_datas_wREZ(mn_data, scenario, argz["scenario_years"], argz["dr"])#sum of years must equal total
    #[println(k*" "*b*" "*string(br["construction_cost"])) for (k,nw) in mn_data["nw"] for (b,br) in nw["ne_branch"]];println()
    mn_data = npvs_costs_datas_4mip(mn_data, scenario, argz["scenario_years"], argz["dr"])#future investment at y scaled to year y=0
    return mn_data
end


#Organizes nw numbers per scenario-year
function multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario = Dict{String, Any}("hours" => ls,"years" => length(scenario_years), "sc_names" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    all_scenario_data=DataFrames.DataFrame()
    #set problem dimension
    dim = scenario["hours"] * length(scenario_years) * length(scenario_names)
    for _nm in scenario_names; push!(scenario["sc_names"], _nm=> Dict{String, Any}());for _yr in scenario_years; push!(scenario["sc_names"][_nm], _yr=> []);end;end

    for (s,(_sc, data_by_sc)) in enumerate(scenario_data);
        data["scenario"][string(s)] = Dict()
        data["scenario_prob"][string(s)] = 1/(length(scenario_names))
        for (t,(_yr, data_by_yr)) in enumerate(data_by_sc);
            all_scenario_data=vcat(all_scenario_data,scenario_data[_sc][_yr])
            start_idx=(s-1)*scenario["hours"]*length(scenario_years)
            start_idx=start_idx+(t-1)*scenario["hours"]
            for h in 1 : scenario["hours"]
                network = start_idx + h
                h2=h+(t-1)*scenario["hours"]
                data["scenario"][string(s)]["$h2"] = network
                push!(scenario["sc_names"][_sc][_yr],network)
            end
        end;
    end
    return all_scenario_data,data,scenario, dim
end


#loads generator cost and profile time series In multi-period simulation
function create_profile_sets_mesh(number_of_hours, data_orig, zs_data, zs, inf_grid, owpp_mva)
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    data=Dict{String,Any}();data["gen"]=Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    data["gen"]=sort!(OrderedCollections.OrderedDict(data_orig["gen"]), by=x->parse(Int64,x))
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end
    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
        for (g, gen) in data["gen"]
            if (gen["type"]>0)#market generator onshore
                extradata["gen"][g]["pmax"][1, d] = inf_grid/pu
                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [(zs_data[!,"EUR_da"*zs[1][gen["gen_bus"]]][d])/e2me,0]
            else#wind gen
                extradata["gen"][g]["pmax"][1, d] = (zs_data[!,"Wnd_MWh"*zs[2][gen["gen_bus"]-length(zs[1])]][d])*owpp_mva[gen["gen_bus"]-length(zs[1])]/pu

                #extradata["gen"][g]["pmax"][1, d]
                #zs_data[!,"Wnd_MWh"*zs[2][gen["gen_bus"]-length(zs[1])]][d]
                #owpp_mva[gen["gen_bus"]-length(zs[1])]/pu

                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [0,0]
            end
        end
    end
    #add loads
    loads=Dict{String,Any}()
    num_of_gens=length(data["gen"])
    for (g, gen) in sort!(OrderedCollections.OrderedDict(extradata["gen"]), by=x->parse(Int64,x))
        if (data["gen"][g]["type"]>0)#market generator onshore
            load=deepcopy(data["gen"][g])
            load["index"]=num_of_gens+1
            load["source_id"][2]=num_of_gens+1
            load["pmin"]=deepcopy(load["pmax"])*-1
            load["pmax"]=0
            push!(loads,string(num_of_gens+1)=>deepcopy(load))
            num_of_gens=num_of_gens+1
        end
    end
    for (l, load) in loads
        extradata["gen"][l] = Dict{String,Any}()
        extradata["gen"][l]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][l]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][l]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
        for (l, load) in loads
            if (load["type"]>0)#market generator onshore
                extradata["gen"][l]["pmax"][1, d] = 0
                extradata["gen"][l]["pmin"][1, d] = (inf_grid/pu)*-1
                extradata["gen"][l]["cost"][d] = [(zs_data[!,"EUR_da"*zs[1][load["gen_bus"]]][d])/e2me,0]
                push!(data_orig["gen"],l=>load)
            else#wind gen
            end
        end
    end

    #set ["type"]
    for (g, gen) in data_orig["gen"]
        gen["type"]=0
    end
    return extradata,data_orig
end

######################## Scaling cost data ###########################
#scale investment to hourly cost spread over the year
function scale_cost_data_2hourly!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / (scenario["hours"])) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / (scenario["hours"])) * x # scale total costs to the planning horizon

    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "invest", rescale_total)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", rescale_total)
    end
end


#scales data to hourly cost spread over horizon years
function scale_cost_data_2yearly!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / (scenario["hours"]*scenario["years"])) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / (scenario["hours"]*scenario["years"])) * x # scale total costs to the planning horizon

    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "invest", rescale_total)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", rescale_total)
    end
end


##################### NPV calculations ###############################
#translates costs (array) to yearly NPV value
function npvs_costs_datas_wREZ(data, scenario, _yrs, _dr)
    _scs=data["scenario"]
    base_yr=parse(Int64,_yrs[1])
    _hrs=deepcopy(scenario["hours"])
    #_yrs=[k for k in keys(scenario["sc_names"][_scs[1]])]
    for (_sci,_sc) in _scs
        _sc_temp=sort!(OrderedCollections.OrderedDict(deepcopy(_sc)), by=x->parse(Int64,x))
        _sc_first=first(_sc_temp)[2]
        _yr=1
        for (_str,_num) in _sc_temp
            if (_num<=_sc_first+_yr*_hrs-1)
                data["nw"][string(_num)]=deepcopy(npv_cost_data_wREZ(deepcopy(data["nw"][string(_num)]),base_yr,parse(Int64,_yrs[_yr]),scenario["planning_horizon"], _dr))
                if (_num==_sc_first+_yr*_hrs-1)
                    _yr=_yr+1
                end
            end
        end
    end
    return data
end


#translates cost (single) to yearly NPV value
function npv_cost_data_wREZ(data,base_yr,current_yr,_ph,_dr::Float64=0.04)
    #println(((_ph-(current_yr-base_yr))/_ph))
    npv_yearly = x -> (1 / (1+_dr)^(current_yr-base_yr)) * x *((_ph-(current_yr-base_yr))/_ph)# npv
    npv_hourly = x -> (1 / (1+_dr)^(current_yr-base_yr)) * x# npv
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        _PM._apply_func!(gen, "cost", npv_hourly)
    end
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        _PM._apply_func!(gen, "invest", npv_yearly)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", npv_yearly)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv_yearly)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", npv_yearly)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", npv_yearly)
        _PM._apply_func!(strg, "inst_cost", npv_yearly)
        _PM._apply_func!(strg, "cost_abs", npv_hourly)
        _PM._apply_func!(strg, "cost_inj", npv_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv_yearly)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv_yearly)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", npv_yearly)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", npv_yearly)
    end
    return data
end


#ensures binary candidates (array) costs sum to proper NPV value over the number of years
function npvs_costs_datas_4mip(data, scenario, _yrs, _dr)
    _scs=data["scenario"]
    _hrs=deepcopy(scenario["hours"])
    for (_sci,_sc) in _scs
        _sc=sort!(OrderedCollections.OrderedDict(_sc), by=x->parse(Int64,x))
        _sc_first=first(_sc)[2]
        for (_str,_num) in _sc
            if (_num<=_sc_first+(length(_yrs)-1)*_hrs-1)
                data["nw"][string(_num)]=deepcopy(npv_cost_data_4mip(deepcopy(data["nw"][string(_num)]),data["nw"][string(_num+_hrs)]))
            end
        end
    end
    return data
end

#ensures binary candidate (single) cost sum to proper NPV value over the number of years
function npv_cost_data_4mip(data0,data1)
    function mip(cost0, cost1)
        return cost0-cost1
    end

    for (b, branch) in get(data0, "ne_branch", Dict{String,Any}())
        data0["ne_branch"][b]["construction_cost"] = mip(data0["ne_branch"][b]["construction_cost"],data1["ne_branch"][b]["construction_cost"])
    end
    for (b, branch) in get(data0, "branchdc_ne", Dict{String,Any}())
        data0["branchdc_ne"][b]["cost"] = mip(data0["branchdc_ne"][b]["cost"],data1["branchdc_ne"][b]["cost"])
    end
    for (c, conv) in get(data0, "convdc_ne", Dict{String,Any}())
        data0["convdc_ne"][c]["cost"] = mip(data0["convdc_ne"][c]["cost"],data1["convdc_ne"][c]["cost"])
    end
    return data0
end

#returns vector of NPV investment limits per year  
function max_invest_per_year(argz)
    max_invest=Float64[]
    for _yr in argz["scenario_years"]
        push!(max_invest,argz["yearly_investment"]/((1+argz["dr"])^(parse(Int64,_yr)-parse(Int64,argz["scenario_years"][1]))));end
    return max_invest
end
######################################################################