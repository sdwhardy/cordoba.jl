#ACDC problem with storage main setup logic
function main_ACDC_wstrg(rt_ex,argz, s)
    ################# Load topology files ###################################
    data, ics_ac, ics_dc, nodes = topology_df(rt_ex, s["relax_problem"], s["AC"])
    ############### defines size and market of genz and wfs ###################
    infinite_grid, genz, wfz, markets = genz_n_wfs(argz["owpp_mva"],nodes,data["baseMVA"])
    push!(argz,"genz"=>genz)
    push!(argz,"wfz"=>wfz)
    #################### Calculates cable options for AC lines
    data=AC_cable_options(data,argz["candidate_ics_ac"],ics_ac,data["baseMVA"])
    print_topology_data_AC(data,markets)#print to verify
    #################### Calculates cable options for DC lines
    data=DC_cable_options(data,argz["candidate_ics_dc"],ics_dc,data["baseMVA"])
    additional_params_PMACDC(data)
    print_topology_data_DC(data,markets)#print to verify
    ##################### load time series data ##############################
    scenario_data, ls = load_time_series(rt_ex,argz)
    #display(scenario_data["EU19"]["2020"])
    push!(argz,"ls"=>ls)
    ##################### multi period setup #################################
    mn_data = multi_period_setup(ls, scenario_data, data, markets, infinite_grid, argz)
    s=update_settings(s, argz, data)
#    if (haskey(s,"home_market") && length(s["home_market"])>0)
#        mn_data=zonal_adjust(mn_data, s);end
    return mn_data, data, argz, s
end

#=
function zonal_adjust(mn_data, s)
    for (k_nw,nw) in mn_data["nw"]
        for (k_c,c_dc) in nw["branchdc_ne"]
            if (issubset([c_dc["fbusdc"]],s["home_market"]) && issubset([c_dc["tbusdc"]],s["home_market"]))
                c_dc["rateA"]=c_dc["rateB"]=c_dc["rateC"]=c_dc["rateA"]*(1-s["balancing_reserve"])
            end
        end
    end
    for (k_nw,nw) in mn_data["nw"]
        for (k_c,c_dc) in nw["ne_branch"]
            if (issubset([c_dc["f_bus"]],s["home_market"]) && issubset([c_dc["t_bus"]],s["home_market"]))
                c_dc["rate_a"]=c_dc["rate_b"]=c_dc["rate_c"]=c_dc["rate_a"]*(1-s["balancing_reserve"])
            end
        end
    end
    return mn_data
end
=#

function update_settings(s, argz, data)
    s["genz"]=argz["genz"]
    s["wfz"]=argz["wfz"]
    s["ic_lim"]=argz["conv_lim"]/data["baseMVA"]
    s["rad_lim"]=maximum([b["rate_a"] for (k,b) in data["ne_branch"]])
    s["scenarios_length"] = length(argz["scenario_names"])
    s["years_length"] = length(argz["scenario_years"])
    s["hours_length"] = argz["ls"]
    return s
end

function max_invest_per_year(argz)
    max_invest=Float64[]
    for _yr in argz["scenario_years"]
        push!(max_invest,argz["yearly_investment"]/((1+argz["dr"])^(parse(Int64,_yr)-parse(Int64,argz["scenario_years"][1]))));end
    return max_invest
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

#adds DC grid to PMACDC
function additional_params_PMACDC(data)
    _PMACDC.process_additional_data!(data)#add extra DC model data
    converter_parameters_rxb(data)#sets converter parameters for loss calc
end

#seperates wfs from genz and defines markets/wfs zones
function genz_n_wfs(owpp_mva,nodes,pu)
    infinite_grid=sum(owpp_mva)*3
    markets_wfs=[String[],String[]]#UK,DE,DK must be in same order as .m file gens
    for (k,cunt) in enumerate(nodes["country"])
        if (nodes["type"][k]>0)
        push!(markets_wfs[1],cunt);else
        push!(markets_wfs[2],cunt);end
    end
    genz=[];wfz=[]
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid/pu));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid/pu));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]/pu));end
    return infinite_grid, genz, wfz, markets_wfs
end

#ensures binary candidates (array) costs sum to proper NPV value over the number of years
function npvs_costs_datas_4mip(data, scenario, _yrs, _dr)
    _scs=data["scenario"]
    _hrs=deepcopy(scenario["hours"])
    for (_sci,_sc) in _scs
        _sc=sort!(OrderedCollections.OrderedDict(_sc),byvalue=true)
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

#translates costs (array) to yearly NPV value
function npvs_costs_datas(data, scenario, _yrs)
    _scs=data["scenario"]
    base_yr=parse(Int64,_yrs[1])
    _hrs=deepcopy(scenario["hours"])
    #_yrs=[k for k in keys(scenario["sc_names"][_scs[1]])]
    for (_sci,_sc) in _scs
        _sc_temp=sort!(OrderedCollections.OrderedDict(deepcopy(_sc)),byvalue=true)
        _sc_first=first(_sc_temp)[2]
        _yr=1
        for (_str,_num) in _sc_temp
            if (_num<=_sc_first+_yr*_hrs-1)
                data["nw"][string(_num)]=deepcopy(npv_cost_data(deepcopy(data["nw"][string(_num)]),base_yr,parse(Int64,_yrs[_yr])))
                if (_num==_sc_first+_yr*_hrs-1)
                    _yr=_yr+1
                end
            end
        end
    end
    return data
end

#translates cost (single) to yearly NPV value
function npv_cost_data(data,base_yr,current_yr,_dr::Float64=0.04)
    npv = x -> (1 / (1+_dr)^(current_yr-base_yr)) * x# npv
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        _PM._apply_func!(gen, "cost", npv)
    end
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        _PM._apply_func!(gen, "invest", npv)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", npv)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", npv)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", npv)
        _PM._apply_func!(strg, "inst_cost", npv)
        _PM._apply_func!(strg, "cost_abs", npv)
        _PM._apply_func!(strg, "cost_inj", npv)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", npv)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", npv)
    end
    return data
end

#translates costs (array) to yearly NPV value
function npvs_costs_datas_wREZ(data, scenario, _yrs, _dr)
    _scs=data["scenario"]
    base_yr=parse(Int64,_yrs[1])
    _hrs=deepcopy(scenario["hours"])
    #_yrs=[k for k in keys(scenario["sc_names"][_scs[1]])]
    for (_sci,_sc) in _scs
        _sc_temp=sort!(OrderedCollections.OrderedDict(deepcopy(_sc)),byvalue=true)
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

#Converts parameters of a converter - left over from Jay Dave
function converter_parameters_rxb(data)
    for (c,conv) in data["convdc_ne"]

        bus = conv["busac_i"]
        base_kV = data["bus"]["$bus"]["base_kv"]
        base_S = sqrt((100*conv["Pacmax"])^2+(100*conv["Qacmax"])^2) #base MVA = 100
        base_Z = base_kV^2/base_S # L-L votlage/3 phase power
        base_Y= 1/base_Z
        conv["xtf"] = 0.10*100/base_S #new X =old X *(100MVA/old Sbase)
        conv["rtf"] = conv["xtf"]/100
        conv["bf"] = 0.08*base_S/100
        conv["xc"] = 0.07*100/base_S #new X =old X *(100MVA/old Zbase)
        conv["rc"] = conv["xc"]/100 #new X =old X *(100MVA/old Zbase)
        rtf = conv["rtf"]
        xtf = conv["xtf"]
        bf = conv["bf"]
        Pmax = conv["Pacmax"]
        Pmin =  conv["Pacmin"]
        Qmax = conv["Qacmax"]
        Qmin =  conv["Qacmin"]

        conv["Imax"] = sqrt(Pmax^2+Qmax^2)
        xc = conv["xc"]
        rc = conv["rc"]
        Imax = conv["Imax"]

    end
    for (c,conv) in data["convdc"]

        bus = conv["busac_i"]
        base_kV = data["bus"]["$bus"]["base_kv"]
        base_S = sqrt((100*conv["Pacmax"])^2+(100*conv["Qacmax"])^2) #base MVA = 100
        base_Z = base_kV^2/base_S # L-L votlage/3 phase power
        base_Y= 1/base_Z
        conv["xtf"] = 0.10*100/base_S #new X =old X *(100MVA/old Sbase)
        conv["rtf"] = conv["xtf"]/100
        conv["bf"] = 0.08*base_S/100
        conv["xc"] = 0.07*100/base_S #new X =old X *(100MVA/old Zbase)
        conv["rc"] = conv["xc"]/100 #new X =old X *(100MVA/old Zbase)
        rtf = conv["rtf"]
        xtf = conv["xtf"]
        bf = conv["bf"]
        Pmax = conv["Pacmax"]
        Pmin =  conv["Pacmin"]
        Qmax = conv["Qacmax"]
        Qmin =  conv["Qacmin"]

        conv["Imax"] = sqrt(Pmax^2+Qmax^2)
        xc = conv["xc"]
        rc = conv["rc"]
        Imax = conv["Imax"]

    end
end

#Sets DC candidate dictionaries with desired candidate qualities
function additional_candidatesICS_DC(data,candidates,ic_data)
    #DC, IC
    ics=[]
    data["branchdc_ne"]=sort!(OrderedCollections.OrderedDict(data["branchdc_ne"]), by=x->parse(Int64,x))
    for (i,dcb) in data["branchdc_ne"]; push!(ics,deepcopy(dcb));end

    data["branchdc_ne"]=Dict{String,Any}()
    for (i,ic) in enumerate(ics)
        for j=1:1:length(candidates);
            ic["source_id"][2]=j+length(candidates)*(i-1);
            ic["index"]=j+length(candidates)*(i-1);
            ic["rateA"]=candidates[j]*first(ic_data[i]);
            ic["length"]=last(ic_data[i])[1];
            push!(data["branchdc_ne"],string(j+length(candidates)*(i-1))=>deepcopy(ic));
        end
    end
    return data
end

#Sets AC candidate dictionaries with desired candidate qualities
function additional_candidatesICS_AC(data,candidates,ic_data)
    #DC, IC
    ics=[]
    data["ne_branch"]=sort!(OrderedCollections.OrderedDict(data["ne_branch"]), by=x->parse(Int64,x))
    for (i,dcb) in data["ne_branch"]; push!(ics,deepcopy(dcb));end

    data["ne_branch"]=Dict{String,Any}()
    for (i,ic) in enumerate(ics)
        for j=1:1:length(candidates);
            ic["source_id"][2]=j+length(candidates)*(i-1);
            ic["index"]=j+length(candidates)*(i-1);
            ic["rate_a"]=candidates[j]*first(ic_data[i]);
            ic["length"]=last(ic_data[i])[1];
            push!(data["ne_branch"],string(j+length(candidates)*(i-1))=>deepcopy(ic));
        end
    end
    return data
end



#for each DC candidate capacity an appropriate cable is selected and characteristics stored
function candidateIC_cost_impedance_DC(bdc,z_base)
    cb=DC_cbl(bdc["rateA"], bdc["length"])
    #bdc["cost"]=cb.costs.cpx_i+cb.costs.cpx_p
    bdc["cost"]=cb.costs.ttl
    #bdc["r"]=((cb.elec.ohm*10^3/cb.num)*cb.length)/z_base
    bdc["r"]=((cb.elec.ohm/cb.num)*cb.length)/z_base
    bdc["rateC"]=bdc["rateB"]=bdc["rateA"]=cb.num*cb.elec.mva
    return bdc
end

function convex2mip_DC(result_mip, data)
    dc_cables=Dict{String,Any}()
    for (j,dc_br) in result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branchdc"]
        dc_min=[];dc_max=[]
        for (k,dc_br_ne) in data["branchdc_ne"]
            if (dc_br_ne["fbusdc"]==data["branchdc"][j]["fbusdc"] && dc_br_ne["tbusdc"]==data["branchdc"][j]["tbusdc"])
                dif=dc_br["p_rateA"]-dc_br_ne["rateA"]

                if (dif>0)
                    push!(dc_min,(k,dif))
                else
                    push!(dc_max,(k,dif))
                end
            end
        end
        if (length(last.(dc_min))>0)
            push!(dc_cables, first.(dc_min)[argmin(last.(dc_min))] => data["branchdc_ne"][first.(dc_min)[argmin(last.(dc_min))]])
        end
        if (length(last.(dc_max))>0)
            push!(dc_cables, first.(dc_max)[argmax(last.(dc_max))] => data["branchdc_ne"][first.(dc_max)[argmax(last.(dc_max))]])
        end
    end
    dc_cables=unique_candidateIC_DC(dc_cables)
    #=dc_cables2=Dict{String,Any}()
    for (i,(k, dc_c)) in enumerate(sort(OrderedCollections.OrderedDict(dc_cables), by=x->parse(Int64,x)))
        dc_c["source_id"][2]=i
        dc_c["status"]=1
        push!(dc_cables2,string(i)=>dc_c)
    end=#
    return dc_cables
end

function convex2mip_AC(result_mip, data)
    ac_cables=Dict{String,Any}()
    for (j,ac_br) in result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branch"]
        ac_min=[];ac_max=[]
        for (k,ac_br_ne) in data["ne_branch"]
            if (ac_br_ne["f_bus"]==data["branch"][j]["f_bus"] && ac_br_ne["t_bus"]==data["branch"][j]["t_bus"])
                dif=ac_br["p_rateAC"]-ac_br_ne["rate_a"]

                if (dif>0)
                    push!(ac_min,(k,dif))
                else
                    push!(ac_max,(k,dif))
                end
            end
        end
        if (length(last.(ac_min))>0)
            push!(ac_cables, first.(ac_min)[argmin(last.(ac_min))] => data["ne_branch"][first.(ac_min)[argmin(last.(ac_min))]])
        end
        if (length(last.(ac_max))>0)
            push!(ac_cables, first.(ac_max)[argmax(last.(ac_max))] => data["ne_branch"][first.(ac_max)[argmax(last.(ac_max))]])
        end
    end
    ac_cables=unique_candidateIC_AC(ac_cables)
    #=ac_cables2=Dict{String,Any}()
    for (i,(k, ac_c)) in enumerate(sort(OrderedCollections.OrderedDict(ac_cables), by=x->parse(Int64,x)))
        ac_c["source_id"][2]=i
        ac_c["br_status"]=1
        push!(ac_cables2,string(i)=>ac_c)
    end=#
    return ac_cables
end

function convex2mip(result_mip, data, mn_data, s)
    s["agent"]=""
    s["relax_problem"]=false
    if (s["AC"]=="1")
        data["ne_branch"]=convex2mip_AC(result_mip, data);end
    data["branchdc_ne"]=convex2mip_DC(result_mip, data)

    for (k,nw) in mn_data["nw"]
        br_nes=deepcopy(data["ne_branch"])
        br_dc_nes=deepcopy(data["branchdc_ne"])
        br_nes2=Dict{String,Any}();br_dc_nes2=Dict{String,Any}()
        for (i,(k,br_ne)) in enumerate(br_nes)
            br_ne["construction_cost"]=nw["ne_branch"][k]["construction_cost"]
            br_ne["source_id"][2]=i
            br_ne["br_status"]= s["AC"]=="0" ? 0 : 1
            push!(br_nes2,string(i)=>br_ne)
        end
        for (i,(k,br_dc_ne)) in enumerate(br_dc_nes)
            br_dc_ne["cost"]=nw["branchdc_ne"][k]["cost"]
            br_dc_ne["source_id"][2]=i
            br_dc_ne["br_status"]=1
            push!(br_dc_nes2,string(i)=>br_dc_ne)
        end
        nw["ne_branch"]=deepcopy(br_nes2)
        nw["branchdc_ne"]=deepcopy(br_dc_nes2)
        for (i,br) in nw["branch"]
            br["br_status"]=0
        end
        for (i,br) in nw["branchdc"]
            br["status"]=0
        end
    end
    data["ne_branch"]=mn_data["nw"]["1"]["ne_branch"]
    data["branchdc_ne"]=mn_data["nw"]["1"]["branchdc_ne"]
    return mn_data, data, s
end
#for each AC candidate capacity an appropriate cable is selected and characteristics stored
function candidateIC_cost_impedance_AC(bac,z_base,s_base)
    cb=AC_cbl(bac["rate_a"], bac["length"])
    bac["construction_cost"]=cb.costs.ttl
    bac["br_r"]=((cb.elec.ohm/cb.num)*cb.length)/z_base
    bac["br_x"]=((cb.elec.xl/cb.num)*cb.length)/z_base
    bac["rate_c"]=bac["rate_b"]=bac["rate_a"]=(cb.num*cb.elec.mva)/s_base
    return bac
end

#ensures that the only candidates considered are unique cable sizes
function unique_candidateIC_DC(cand_ics)
    copy_cand_ics=deepcopy(cand_ics)
    for (i,dcb) in cand_ics
        for (j,tdcb) in copy_cand_ics
            if (i!=j && dcb["fbusdc"]==tdcb["fbusdc"] && dcb["tbusdc"]==tdcb["tbusdc"] &&  isapprox(dcb["rateA"],tdcb["rateA"]; atol = 1))
                delete!(copy_cand_ics,j)
                break
            end
        end
    end
    return copy_cand_ics
end

function unique_candidateIC_AC(cand_ics)
    copy_cand_ics=deepcopy(cand_ics)
    for (i,dcb) in cand_ics
        for (j,tdcb) in copy_cand_ics
            if ((i!=j && dcb["f_bus"]==tdcb["f_bus"] && dcb["t_bus"]==tdcb["t_bus"] &&  isapprox(dcb["rate_a"],tdcb["rate_a"]; atol = 1)) || isapprox(tdcb["rate_a"],0; atol = 1))
                delete!(copy_cand_ics,j)
                break
            end
        end
    end
    return copy_cand_ics
end

function cable_reorder_rename(cbls)
    cbls2=sort(OrderedCollections.OrderedDict(deepcopy(cbls)), by = x -> parse(Int64,x))
    ks_orig=keys(cbls2)
    for (k,k_orig) in enumerate(ks_orig)
        if (string(k_orig)!=string(k))
            cbls[string(k)]=cbls[k_orig];
            delete!(cbls,string(k_orig))
        end
    end
    return cbls
end

#divides time series into 24 hour groups
function daily_tss(ts)
    ts_daily=ts[1:24]
    for i=25:24:length(ts)-24
        ts_daily=hcat(ts_daily,ts[i:i+23])
    end
    return ts_daily
end

#divides time series in n hour groups
function half_daily_tss(ts,n)
    ts_daily=ts[1:n]
    for i=(n+1):n:length(ts)-n
        ts_daily=hcat(ts_daily,ts[i:i+(n-1)])
    end
    return ts_daily
end

#reads time series from file into dataframe for a given set of scenarios and years
function get_scenario_year_tss(sc_nms,sc_yrs)
    scenario_data = Dict{String,Any}()
    for _sc in sc_nms
        push!(scenario_data,_sc=>Dict{String,Any}())
        for _yr in sc_yrs
            df=CSV.read("./test/data/input/scenarios/convex_wBE/"*_sc*_yr*".csv", DataFrames.DataFrame)
            push!(scenario_data[_sc],_yr=>df)
        end
    end
    return scenario_data
end

#Organizes nw numbers per scenario-year
function multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario = Dict{String, Any}("hours" => ls,"years" => length(scenario_years), "sc_names" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    all_scenario_data=DataFrame()
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

#used to run a sim wjile fixing the values of candidate cables and converters
function fix_cables_and_converters(mn_data_nw, fixed_variables,cables, converters)
    for (key,nw) in mn_data_nw
        push!(fixed_variables,key=>Dict{String,Any}())
        push!(fixed_variables[key],"baseMVA" => nw["baseMVA"])
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"branchdc_ne" => Dict{String, Any}())
        for (key_br,br) in nw["branchdc_ne"]
            push!(fixed_variables[key]["branchdc_ne"],key_br => Dict{String, Any}())
            if (issubset([parse(Int64,key_br)],cables))
                push!(fixed_variables[key]["branchdc_ne"][key_br],"isbuilt" => 1.00)
            else
                push!(fixed_variables[key]["branchdc_ne"][key_br],"isbuilt" => 0.00)
            end
        end

        #if needed add convdc_ne here
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"convdc" => Dict{String, Any}())
        for (key_c,c) in nw["convdc"]
            push!(fixed_variables[key]["convdc"],key_c => Dict{String, Any}())
            push!(fixed_variables[key]["convdc"][key_c],"p_pacmax" => converters[parse(Int64,key_c)])
        end
    end
    return fixed_variables
end


###################################################### DEPRICATED ##########################################
#############################################
#storage setup for Ancillary services analysis
#not verified and will likely be depricated as WP2 is handling this
function add_storage_profile(dim, data, extradata, zs_data, zs, number_of_hours)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU

    storage=[(i,b) for (i,b) in data["ne_storage"]]
    sort!(storage, by=x->x[2]["energy_rating"])
    extradata["ne_storage"] = Dict{String,Any}()
    for (b, bat) in data["ne_storage"]
        extradata["ne_storage"][b] = Dict{String,Any}()
        extradata["ne_storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["charge_rating"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["discharge_rating"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        #up_reg

        up_pe=[(zs_data["EUR_up"*z][d],zs_data["MWh_up"*z][d]) for z in zs[2]]
        sort!(up_pe, by = x -> x[1], rev=true)
        up_pe=(first.(up_pe)./e2me,last.(up_pe)./pu)

        dwn_pe=[(zs_data["EUR_dwn"*z][d],zs_data["MWh_dwn"*z][d]) for z in zs[2]]
        elec_prices=[(zs_data["EUR_id"*z][d],Inf) for z in zs[2]]
        sort!(elec_prices, by = x -> x[1])
        #push!(dwn_pe,elec_prices[1])
        sort!(dwn_pe, by = x -> x[1])
        dwn_pe=(first.(dwn_pe)./e2me,last.(dwn_pe)./pu)


        for (b,bat) in storage
            #Set discharge rating and cost
            dcr_temp=deepcopy(bat["discharge_rating"])
            extradata["ne_storage"][b]["discharge_rating"][1, d] = 0
            extradata["ne_storage"][b]["cost_inj"][1, d] = 0
            for upe=1:1:length(up_pe[1])
                if (dcr_temp<=last(up_pe)[upe])
                    extradata["ne_storage"][b]["discharge_rating"][1, d] = extradata["ne_storage"][b]["discharge_rating"][1, d] + dcr_temp
                    extradata["ne_storage"][b]["cost_inj"][1, d] = extradata["ne_storage"][b]["cost_inj"][1, d] + first(up_pe)[upe]
                    last(up_pe)[upe]=last(up_pe)[upe]-dcr_temp
                    break
                elseif (last(up_pe)[upe]>1e-3 && last(up_pe)[upe]<dcr_temp)
                    extradata["ne_storage"][b]["discharge_rating"][1, d] = extradata["ne_storage"][b]["discharge_rating"][1, d] + last(up_pe)[upe]
                    extradata["ne_storage"][b]["cost_inj"][1, d] = ((extradata["ne_storage"][b]["discharge_rating"][1, d]-last(up_pe)[upe])/extradata["ne_storage"][b]["discharge_rating"][1, d])*extradata["ne_storage"][b]["cost_inj"][1, d] + (last(up_pe)[upe]/extradata["ne_storage"][b]["discharge_rating"][1, d])*first(up_pe)[upe]
                    dcr_temp=dcr_temp-last(up_pe)[upe]
                    last(up_pe)[upe] = 0
                elseif (last(up_pe)[upe]<=1e-3 && extradata["ne_storage"][b]["discharge_rating"][1, d] <= 1e-3)
                    extradata["ne_storage"][b]["discharge_rating"][1, d] = deepcopy(bat["discharge_rating"])
                    extradata["ne_storage"][b]["cost_inj"][1, d] = 0
                    break
                end
            end
            #Set charge rating and cost
            cr_temp=deepcopy(bat["charge_rating"])
            extradata["ne_storage"][b]["charge_rating"][1, d] = 0
            extradata["ne_storage"][b]["cost_abs"][1, d] = 0
            for upe=1:1:length(dwn_pe[1])
                if (cr_temp<=last(dwn_pe)[upe] && last(dwn_pe)[upe]!=Inf)
                    extradata["ne_storage"][b]["charge_rating"][1, d] = extradata["ne_storage"][b]["charge_rating"][1, d] + cr_temp
                    extradata["ne_storage"][b]["cost_abs"][1, d] = extradata["ne_storage"][b]["cost_abs"][1, d] + first(dwn_pe)[upe]
                    last(dwn_pe)[upe]=last(dwn_pe)[upe]-cr_temp
                    break
                elseif (last(dwn_pe)[upe]>1e-3 && last(dwn_pe)[upe]<cr_temp)
                    extradata["ne_storage"][b]["charge_rating"][1, d] = extradata["ne_storage"][b]["charge_rating"][1, d] + last(dwn_pe)[upe]
                    extradata["ne_storage"][b]["cost_abs"][1, d] = ((extradata["ne_storage"][b]["charge_rating"][1, d]-last(dwn_pe)[upe])/extradata["ne_storage"][b]["charge_rating"][1, d])*extradata["ne_storage"][b]["cost_abs"][1, d] + (last(dwn_pe)[upe]/extradata["ne_storage"][b]["charge_rating"][1, d])*first(dwn_pe)[upe]
                    cr_temp=cr_temp-last(dwn_pe)[upe]
                    last(dwn_pe)[upe] = 0
                elseif ((last(dwn_pe)[upe]<=1e-3 || last(dwn_pe)[upe]==Inf) && extradata["ne_storage"][b]["charge_rating"][1, d] <= 1e-3)
                    extradata["ne_storage"][b]["charge_rating"][1, d] = deepcopy(bat["charge_rating"])
                    extradata["ne_storage"][b]["cost_abs"][1, d] = first(dwn_pe)[upe]
                    break
                end
            end
        end
    end

    return extradata,data
end

#
#=
function candidateIC_cost(bdc)
    println("from: "*string(bdc["fbusdc"])*" to: "*string(bdc["tbusdc"]))
    bdc["cost"],bdc["r"],bdc["rateA"],bdc["r"]=dc_cable_cost_impedance(bdc["rateA"],bdc["length"])
    #=if (bdc["rateC"]==-1)#on-on
        println("from: "*string(bdc["fbusdc"])*" to: "*string(bdc["tbusdc"]))
        bdc["cost"],bdc["r"],bdc["rateA"]=on_on_ic(bdc["rateA"],bdc["rateB"])
    elseif (bdc["rateC"]==0)#off-off
        println("from: "*string(bdc["fbusdc"])*" to: "*string(bdc["tbusdc"]))
        bdc["cost"],bdc["r"],bdc["rateA"]=off_off_ic(bdc["rateA"],bdc["rateB"])
    elseif (bdc["rateC"]==1)#on-off
        println("from: "*string(bdc["fbusdc"])*" to: "*string(bdc["tbusdc"]))
        bdc["cost"],bdc["r"],bdc["rateA"]=on_off_ic(bdc["rateA"],bdc["rateB"])
    end=#
    bdc["rateC"]=bdc["rateB"]=bdc["rateA"]
    return bdc
end=#

#=
function create_profile_sets_owpps(number_of_hours, data, zs_data, zs, inf_grid, owpp_mva)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
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
                extradata["gen"][g]["cost"][d] = [(zs_data["EUR_da"*zs[gen["gen_bus"]]][d])/e2me,0]
            else#wind gen
                extradata["gen"][g]["pmax"][1, d] = owpp_mva/pu
                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [0,0]
            end
        end
    end
    #add loads
    loads=Dict{String,Any}()
    num_of_gens=length(data["gen"])
    for (g, gen) in extradata["gen"]
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
                extradata["gen"][l]["cost"][d] = [(zs_data["EUR_da"*zs[load["gen_bus"]]][d])/e2me,0]
                push!(data["gen"],l=>load)
            else#wind gen
            end
        end
    end

    #set ["type"]
    for (g, gen) in data["gen"]
        gen["type"]=0
    end
    return extradata,data
end
=#
#=
function dc_cable_cost_impedance(mva,km)
    cb=DC_cbl(mva, km)
    cost=cb.costs.cpx_i+cb.costs.cpx_p
    println("DC connection cost: total "*string(cost)*" mva "*string(cb.num*cb.elec.mva)*" km "*string(km))
    return cost,(cb.elec.ohm/cb.num)*km,cb.num*cb.elec.mva,cb.elec.ohm
end
=#
#
#=function get_scenario_tss(sc_nms,sc_yrs)
    scenario_data = Dict{String,Any}()
    for _yr in sc_yrs
        push!(scenario_data,_yr=>Dict{String,Any}())
        for _sc in sc_nms
            df=CSV.read("./test/data/input/scenarios/scenario_"*_sc*_yr*".csv", DataFrames.DataFrame)
            push!(scenario_data[_yr],_sc=>df)
        end
    end
    return scenario_data
end=#

#=function get_scenario_year_tss(sc_nms,sc_yrs)
    scenario_data = Dict{String,Any}()
    for _sc in sc_nms
        push!(scenario_data,_sc=>Dict{String,Any}())
        for _yr in sc_yrs
            df=CSV.read("./test/data/input/scenarios/scenario_"*_sc*_yr*".csv", DataFrames.DataFrame)
            push!(scenario_data[_sc],_yr=>df)
        end
    end
    return scenario_data
end=#

#=function multi_period_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario = Dict{String, Any}("hours" => ls, "sc_years" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    all_scenario_data=DataFrame()
    #set problem dimension
    scenario["hours"]= ls
    dim = scenario["hours"] * length(scenario_years) * length(scenario_names)
    scene_count=1;for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
        scenario["sc_years"][string(scene_count)] = Dict{String, Any}()
        scenario["sc_years"][string(scene_count)]["year"] = parse(Int64,_yr)#year of data
        scenario["sc_years"][string(scene_count)]["probability"] = 1/(length(scenario_years) * length(scenario_names))
        all_scenario_data=vcat(all_scenario_data,scenario_data[_yr][_sc])
        data["scenario"][string(scene_count)] = Dict()
        data["scenario_prob"][string(scene_count)] = scenario["sc_years"][string(scene_count)]["probability"]
        start_idx=(scene_count-1)*scenario["hours"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][string(scene_count)]["$h"] = network
        end
        scene_count=scene_count+1
    end;end
    return all_scenario_data,data,scenario, dim
end=#

#=
function scale_cost_data_cordoba_convexafy!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / scenario["hours"]) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / scenario["hours"]) * x # scale total costs to the planning horizon
    for (b, branch) in data["branchdc"]
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in data["convdc"]
        _PM._apply_func!(conv, "cost", rescale_total)
    end
end=#
#=
function multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario = Dict{String, Any}("hours" => ls, "sc_years" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    all_scenario_data=DataFrame()
    #set problem dimension
    scenario["hours"]= ls
    dim = scenario["hours"] * length(scenario_years) * length(scenario_names)
    scene_count=1;for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
        scenario["sc_years"][string(scene_count)] = Dict{String, Any}()
        scenario["sc_years"][string(scene_count)]["year"] = parse(Int64,_yr)#year of data
        scenario["sc_years"][string(scene_count)]["probability"] = 1/(length(scenario_years) * length(scenario_names))
        all_scenario_data=vcat(all_scenario_data,scenario_data[_yr][_sc])
        data["scenario"][string(scene_count)] = Dict()
        data["scenario_prob"][string(scene_count)] = scenario["sc_years"][string(scene_count)]["probability"]
        start_idx=(scene_count-1)*scenario["hours"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][string(scene_count)]["$h"] = network
        end
        scene_count=scene_count+1
    end;end
    return all_scenario_data,data,scenario, dim
end=#


#=
function combine_profile_data_sets_entso_scenario(zs,data, n,_sc,_yr, scenario)
    data, zs_data = get_profile_data_sets_entso_scenario(zs,data, n,_sc[1],_yr, scenario)
    unique!(zs_data,:time_stamp)
    data, zs_dataB = get_profile_data_sets_entso_scenario(zs,data, n,_sc[2],_yr, scenario)
    unique!(zs_dataB,:time_stamp)
    data, zs_dataC = get_profile_data_sets_entso_scenario(zs,data, n,_sc[3],_yr, scenario)
    unique!(zs_dataC,:time_stamp)
    filter!(row -> row.time_stamp in zs_data.time_stamp, zs_dataB)
    filter!(row -> row.time_stamp in zs_data.time_stamp, zs_dataC)
    filter!(row -> row.time_stamp in zs_dataB.time_stamp, zs_data)
    filter!(row -> row.time_stamp in zs_dataB.time_stamp, zs_dataC)
    filter!(row -> row.time_stamp in zs_dataC.time_stamp, zs_data)
    filter!(row -> row.time_stamp in zs_dataC.time_stamp, zs_dataB)
    zs_data["EUR_daUK"]=zs_data["EUR_daUK"].*(5/25).+zs_dataB["EUR_daUK"].*(10/25).+zs_dataC["EUR_daUK"].*(10/25)
    zs_data["EUR_daDK"]=zs_data["EUR_daDK"].*(5/25).+zs_dataB["EUR_daDK"].*(10/25).+zs_dataC["EUR_daDK"].*(10/25)
    zs_data["EUR_daDE"]=zs_data["EUR_daDE"].*(5/25).+zs_dataB["EUR_daDE"].*(10/25).+zs_dataC["EUR_daDE"].*(10/25)
    return data,zs_data
end=#

#=function scale_bat_data_cordoba!(data, scenario)
    rescale_hourly = x -> (scenario["hours"] / (8760*scenario["planning_horizon"])) * x # yearly limit on energy absoption
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "max_energy_absorption", rescale_hourly)
    end
end=#

#=for (k0,ss) in scenario_data; for (k1,s) in ss
    s.time_stamp=format_datetime(s.time_stamp);end;end
for (k0,ss) in scenario_data; for (k1,s) in ss
CSV.write("./test/data/input/scenarios/convex/"*string(k0)*string(k1)*".csv",s);end;end
function format_datetime(ts)
    formated_date=[]
    for t in ts
        push!(formated_date,DateTime(string(t[1:16]), dateformat"dd.mm.yyyy HH:MM"))#01.01.2018 00:00
    end
    return formated_date
end=#

#=
function get_profile_data_sets_mesh(zs,data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    zs_data=[]
    for z in zs
        df=CSV.read("./test/data/input/"*z*"data.csv", DataFrames.DataFrame)
        #df=CSV.read("./test/data/input/scenarios/"*yr*"/"*z*"_"*scen*".csv", DataFrames.DataFrame)
        colnames = ["time_stamp","Wnd_MWh"*z,"EUR_da"*z,"EUR_id"*z,"MWh_up"*z,"EUR_up"*z,"MWh_dwn"*z,"EUR_dwn"*z]
        #colnames = ["time_stamp","Wnd_MWh"*z,"EUR_da"*z]
        names!(df, Symbol.(colnames))
        push!(zs_data,df)
    end
    zsd=zs_data[1];for z in zs_data[2:end];
    zsd=innerjoin(zsd,z, makeunique=true,on=:time_stamp);end
    Ytr=PCA_cluster(zsd)
    Ytr = convert(Matrix, zsd[1:1:end,2:end])'
    cluster=kmeans_cluster(Ytr)
    cluster=kmedoid_cluster(Ytr)
    n_samples=n_samps(cluster,n)
    zsd=zsd[n_samples,:];
    s="1";scenario["sc_years"][s]
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]
        scenario["hours"]=length(zsd[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, zsd
end=#

#=function get_profile_data_sets_entso_scenario(zs,data, n,scen,yr, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    zs_data=[]
    for z in zs
        #df=CSV.read("./test/data/input/"*z*"data.csv", DataFrames.DataFrame)
        df=CSV.read("./test/data/input/scenarios/"*yr*"/"*z*"_"*scen*".csv", DataFrames.DataFrame)
        if (haskey(df,"Column1"))
            df=select!(df, Not(:Column1))
        end
        if (haskey(df,"Day-ahead Price [GBP/MWh]"))
            df=select!(df, Not("Day-ahead Price [GBP/MWh]"))
        end
        #colnames = ["time_stamp","Wnd_MWh"*z,"EUR_da"*z,"EUR_id"*z,"MWh_up"*z,"EUR_up"*z,"MWh_dwn"*z,"EUR_dwn"*z]
        colnames = ["time_stamp","Wnd_MWh"*z,"EUR_da"*z]
        names!(df, Symbol.(colnames))
        push!(zs_data,df)
    end
    zsd=zs_data[1];for z in zs_data[2:end];
    zsd=innerjoin(zsd,z, makeunique=true,on=:time_stamp);end
    #=Ytr=PCA_cluster(zsd)
    Ytr = convert(Matrix, zsd[1:1:end,2:end])'
    cluster=kmeans_cluster(Ytr)
    cluster=kmedoid_cluster(Ytr)
    n_samples=n_samps(cluster,n)
    zsd=zsd[n_samples,:];=#
    zsd = zsd[completecases(zsd), :]
    disallowmissing!(zsd)
    s="1";scenario["sc_years"][s]
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]
        scenario["hours"]=length(zsd[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info

    return data, zsd
end=#

#=
function create_profile_sets(number_of_hours, data, df0, df1,ic_mva,owpp_mva)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU
    #e2me=1
    da=0.835;id=0.165
    #e2me=1#into ME/PU
    extradata = Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df0.wind[d]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [da*df0.daprice[d]/e2me+id*df0.idprice[d]/e2me,0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df0.wind[d]
            extradata["gen"]["2"]["cost"][d] = [da*df0.daprice[d]/e2me+id*df0.idprice[d]/e2me,0]

            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df0.wind[d]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [da*df1.daprice[d]/e2me+id*df1.idprice[d]/e2me,0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df0.wind[d]
                extradata["gen"]["4"]["cost"][d] = [da*df1.daprice[d]/e2me+id*df1.idprice[d]/e2me,0]
            #Wind generator
                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/pu*df0.wind[d]
                extradata["gen"]["5"]["pmin"][1, d] = 0
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]
    end
    return extradata
end
=#
#=
function create_profile_sets_wstrg(number_of_hours, data, df0, df1,ic_mva,owpp_mva)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU
    #e2me=1
    da=0.835;id=0.165
    #e2me=1#into ME/PU
    extradata = Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    storage=[(i,b) for (i,b) in data["ne_storage"]]
    sort!(storage, by=x->x[2]["energy_rating"])
    extradata["ne_storage"] = Dict{String,Any}()
    for (b, bat) in data["ne_storage"]
        extradata["ne_storage"][b] = Dict{String,Any}()
        extradata["ne_storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["charge_rating"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["discharge_rating"] = Array{Float64,2}(undef, 1, number_of_hours)
    end

    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df0.wind[d]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [da*df0.daprice[d]/e2me+id*df0.idprice[d]/e2me,0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df0.wind[d]
            extradata["gen"]["2"]["cost"][d] = [da*df0.daprice[d]/e2me+id*df0.idprice[d]/e2me,0]

            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df0.wind[d]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [da*df1.daprice[d]/e2me+id*df1.idprice[d]/e2me,0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df0.wind[d]
                extradata["gen"]["4"]["cost"][d] = [da*df1.daprice[d]/e2me+id*df1.idprice[d]/e2me,0]
            #Wind generator
                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/pu*df0.wind[d]
                extradata["gen"]["5"]["pmin"][1, d] = 0
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]

                #up_reg
                up_prices=[df0.regup_price[d],df1.regup_price[d]];
                best_up_price=deepcopy(findmax(up_prices));up_prices[best_up_price[2]]=Inf
                worst_up_price=findmin(up_prices)
                best_up_mwh=[df0.regup_mwh[d],df1.regup_mwh[d]][best_up_price[2]]
                worst_up_mwh=[df0.regup_mwh[d],df1.regup_mwh[d]][worst_up_price[2]]
                #down_reg
                dwn_prices=[df0.regdwn_price[d],df1.regdwn_price[d]];
                best_dwn_price=deepcopy(findmax(dwn_prices));dwn_prices[best_dwn_price[2]]=Inf
                worst_dwn_price=findmin(dwn_prices)
                best_dwn_mwh=[df0.regdwn_mwh[d],df1.regdwn_mwh[d]][best_dwn_price[2]]
                worst_dwn_mwh=[df0.regdwn_mwh[d],df1.regdwn_mwh[d]][worst_dwn_price[2]]

                best_up_mwh = best_up_mwh>0 ? best_up_mwh/pu : 0
                best_dwn_mwh = best_dwn_mwh>0 ? best_dwn_mwh/pu : 0
                worst_up_mwh = worst_up_mwh>0 ? worst_up_mwh/pu : 0
                worst_dwn_mwh = worst_dwn_mwh>0 ? worst_dwn_mwh/pu : 0

                best_up_price = best_up_price[1]>0 ? best_up_price[1]/e2me : 0
                best_dwn_price = best_dwn_price[1]>0 ? best_dwn_price[1]/e2me : 0
                worst_up_price = worst_up_price[1]>0 ? worst_up_price[1]/e2me : 0
                worst_dwn_price = worst_dwn_price[1]>0 ? worst_dwn_price[1]/e2me : 0

                for (b,bat) in storage
                    dcr_temp=deepcopy(data["ne_storage"][b]["discharge_rating"])
                    bup_mwh_temp=deepcopy(best_up_mwh)
                    wdwn_mwh_temp=deepcopy(worst_up_mwh)

                    best_up_mwh = dcr_temp>best_up_mwh ?  best_up_mwh : dcr_temp#can it be from the battery?
                    best_up_mwh = best_up_mwh>0 ?  best_up_mwh : 0
                    worst_up_mwh = dcr_temp>(best_up_mwh+worst_up_mwh) ?  worst_up_mwh : dcr_temp-best_up_mwh#can additional be from the battery?
                    worst_up_mwh = worst_up_mwh>0 ?  worst_up_mwh : 0
                    if ((best_up_mwh+worst_up_mwh)>0.0001)
                        extradata["ne_storage"][b]["discharge_rating"][1, d]=best_up_mwh+worst_up_mwh
                        extradata["ne_storage"][b]["cost_inj"][1, d] = (best_up_mwh/(best_up_mwh+worst_up_mwh))*best_up_price+((worst_up_mwh)/(best_up_mwh+worst_up_mwh))*worst_up_price
                        best_up_mwh=bup_mwh_temp-best_up_mwh
                        worst_up_mwh=wdwn_mwh_temp-worst_up_mwh
                    else
                        extradata["ne_storage"][b]["discharge_rating"][1, d]=0
                        extradata["ne_storage"][b]["cost_inj"][1, d]=0
                    end
                    if (isnan(extradata["ne_storage"][b]["discharge_rating"][1, d]))
                        println("Bad data - detected at "*string(d)*" - "*string(b))
                        extradata["ne_storage"][b]["discharge_rating"][1, d]=0
                        extradata["ne_storage"][b]["cost_inj"][1, d]=0
                    end
                    cr_temp=deepcopy(data["ne_storage"][b]["charge_rating"])
                    bdwn_mwh_temp=deepcopy(best_dwn_mwh)
                    wdwn_mwh_temp=deepcopy(worst_dwn_mwh)

                    best_dwn_mwh = cr_temp>best_dwn_mwh ?  best_dwn_mwh : cr_temp#can it be from the battery?
                    best_dwn_mwh = best_dwn_mwh>0 ?  best_dwn_mwh : 0
                    worst_dwn_mwh = cr_temp>(best_dwn_mwh+worst_dwn_mwh) ?  worst_dwn_mwh : cr_temp-best_dwn_mwh#can additional be from the battery?
                    worst_dwn_mwh = worst_dwn_mwh>0 ?  worst_dwn_mwh : 0
                    if ((best_dwn_mwh+worst_dwn_mwh)>0.0001)
                        extradata["ne_storage"][b]["charge_rating"][1, d]=cr_temp
                        extradata["ne_storage"][b]["cost_abs"][1, d] = (best_dwn_mwh/cr_temp)*best_dwn_price+(worst_dwn_mwh/cr_temp)*worst_dwn_price
                        best_dwn_mwh=bdwn_mwh_temp-best_dwn_mwh
                        worst_dwn_mwh=wdwn_mwh_temp-worst_dwn_mwh
                    else
                        extradata["ne_storage"][b]["charge_rating"][1, d]=cr_temp
                        extradata["ne_storage"][b]["cost_abs"][1, d]=0
                    end
                end
    end
    return extradata
end
=#
#=
function get_profile_data_sets(d1,d2,data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];

    df0=DataFrame("time_stamp"=>[],"daprice"=>[],"idprice"=>[],"wind"=>[],"regup_price"=>[],"regup_mwh"=>[],"regdwn_price"=>[],"regdwn_mwh"=>[])
    df1=DataFrame("time_stamp"=>[],"daprice"=>[],"idprice"=>[],"wind"=>[],"regup_price"=>[],"regup_mwh"=>[],"regdwn_price"=>[],"regdwn_mwh"=>[])#MWh_up,EUR_up,MWh_dwn,EUR_dwn
    z0_data=CSV.read("./test/data/input/"*d1*".csv", DataFrames.DataFrame)
    z1_data=CSV.read("./test/data/input/"*d2*".csv", DataFrames.DataFrame)
    z01_data=innerjoin(z0_data,z1_data, makeunique=true,on=:time_stamp)
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]

        #tss=[DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.be_eumwh.+0.165.*ukbe_data.be_costid),DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.uk_eumwh.+0.165.*ukbe_data.uk_costid)]
        tss=[DataFrame("time_stamp"=>z01_data.time_stamp, "price"=>abs.((0.835.*z01_data.EUR_da.+0.165.*z01_data.EUR_id).-(0.835.*z01_data.EUR_da_1.+0.165.*z01_data.EUR_id_1))),DataFrame("time_stamp"=>z01_data.time_stamp, "price"=>z01_data.Wnd_MWh)]
        tss_bins=cluster_ts(tss,n)
        sc=sample_cluster(tss_bins,tss,n)
        sort!(sc)
        for t in sc
            #Zone 0
            push!(df0,[t,z0_data[!,:EUR_da][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:EUR_id][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:Wnd_MWh][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:MWh_up][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:EUR_up][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:MWh_dwn][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:EUR_dwn][findfirst(isequal(t),z0_data[!,:time_stamp])]])
            #Zone 1
            push!(df1,[t,z1_data[!,:EUR_da][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:EUR_id][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:Wnd_MWh][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:MWh_up][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:EUR_up][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:MWh_dwn][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:EUR_dwn][findfirst(isequal(t),z1_data[!,:time_stamp])]])
        end

        scenario["hours"]=length(df0[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, df0,df1
end
=#
#=
function get_n_profile_data(data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    if haskey(scenario, "mc")
        monte_carlo = scenario["mc"]
    else
        monte_carlo = false
    end
    df=DataFrame("time_stamp"=>[],"z0_daprice"=>[],"z0_idprice"=>[],"z1_daprice"=>[],"z1_idprice"=>[],"z0_wind"=>[],"z1_wind"=>[])
    ukbe_data=CSV.read("./test/data/cordoba/input/ukbe_ts_2.csv", DataFrames.DataFrame)
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]

        #tss=[DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.be_eumwh.+0.165.*ukbe_data.be_costid),DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.uk_eumwh.+0.165.*ukbe_data.uk_costid)]
        tss=[DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>abs.((0.835.*ukbe_data.be_eumwh.+0.165.*ukbe_data.be_costid).-(0.835.*ukbe_data.uk_eumwh.+0.165.*ukbe_data.uk_costid))),DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>ukbe_data.be_wind)]
        tss_bins=cluster_ts(tss,n)
        sc=sample_cluster(tss_bins,tss,n)
        sort!(sc)
        for t in sc
            push!(df,[t,ukbe_data[!,:be_eumwh][findfirst(isequal(t),ukbe_data[!,:time_stamp])], ukbe_data[!,:be_costid][findfirst(isequal(t),ukbe_data[!,:time_stamp])],ukbe_data[!,:uk_eumwh][findfirst(isequal(t),ukbe_data[!,:time_stamp])], ukbe_data[!,:uk_costid][findfirst(isequal(t),ukbe_data[!,:time_stamp])],ukbe_data[!,:be_wind][findfirst(isequal(t),ukbe_data[!,:time_stamp])], ukbe_data[!,:uk_wind][findfirst(isequal(t),ukbe_data[!,:time_stamp])]])
        end

        scenario["hours"]=length(df[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, df
end
=#



#=function scale_cost_data_cordoba!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / scenario["hours"]) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / scenario["hours"]) * x # scale total costs to the planning horizon
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
        _PM._apply_func!(branch, "co2_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
        _PM._apply_func!(branch, "co2_cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
        _PM._apply_func!(conv, "co2_cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "co2_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end

end=#
