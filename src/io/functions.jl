################################ zonal/nodal market models main function #####################################
function zonal_market_main(s)
    hm=deepcopy(s["home_market"]);
    mn_data, data, s = data_setup_zonal(s);#Build data structure for given options    
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1);#select solver
    result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s);#Solve problem
    #print_solution_wcost_data(result_mip, s, data);
    mn_data, data, s = data_setup_nodal(s);#Build data structure for given options
    mn_data, s = set_inter_zonal_grid(result_mip,mn_data,s);
    s["home_market"]=[]    
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
    result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    #print_solution_wcost_data(result_mip, s, data)
    s["home_market"]=hm
    s["rebalancing"]=true
    s["relax_problem"]=true
    s["output"]["duals"]=true
    mn_data, data, s = data_update(s,result_mip);#Build data structure for given options
    mn_data, s = set_rebalancing_grid(result_mip,mn_data,s);
    s, mn_data= remove_integers(result_mip,mn_data,data,s);
    result_mip_hm_prices = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    s["home_market"]=[]
    mn_data, data, s = data_update(s,result_mip);#Build data structure for given options
    mn_data, s = set_rebalancing_grid(result_mip,mn_data,s);
    s, mn_data= remove_integers(result_mip,mn_data,data,s);
    result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    result_mip= hm_market_prices(result_mip, result_mip_hm_prices)
    return result_mip, data, mn_data, s
end

function nodal_market_main(s)
    mn_data, data, s = data_setup_nodal(s);
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
    result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    print_solution_wcost_data(result_mip, s, data)
    s["rebalancing"]=true
    s["relax_problem"]=true
    s["output"]["duals"]=true
    mn_data, data, s = data_update(s,result_mip);#Build data structure for given options
    mn_data, s = set_rebalancing_grid(result_mip,mn_data,s);
    s, mn_data= remove_integers(result_mip,mn_data,data,s);
    result_mip =  cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem=#
    return result_mip, data, mn_data, s
end
########################## 

##################### Topology input data ############################
#ACDC problem with storage main logic
function main_ACDC_wstrg(rt_ex,argz, s)
    ################# Load topology files ###################################
    topology_df(rt_ex, s["relax_problem"], s["AC"])#creates .m file
    data, ics_ac, ics_dc, nodes = filter_mfile_cables(rt_ex)#loads resulting topology and filters for candidate cables
    ############### defines size and market of genz and wfs ###################
    infinite_grid, genz, wfz, markets = genz_n_wfs(argz["owpp_mva"],nodes,data["baseMVA"],s)
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
    push!(argz,"ls"=>ls)
    ##################### multi period setup #################################
	s=update_settings(s, argz, data)
    mn_data, xd  = multi_period_setup(ls, scenario_data, data, markets, infinite_grid, argz, s)
	s["xd"]=xd
	push!(s,"max_invest_per_year"=>max_invest_per_year(argz))
    return  mn_data, data, argz, s
    #return  scenario_data, data, argz, s, ls, markets, infinite_grid
end

function get_topology_data(s)
    ################# Load topology files ###################################
    topology_df(s["rt_ex"], s["relax_problem"], s["AC"])#creates .m file
    data, s["ics_ac"], s["ics_dc"], s["nodes"] = filter_mfile_cables(s["rt_ex"])#loads resulting topology and filters for candidate cables
    return data, s
end

function get_scenario_data(s)
    ############### defines size and market of genz and wfs ###################
	scenario_data=FileIO.load(s["scenario_data_file"])
    ######## Batteries are removed and modeled seperately time series #########
    for (k_sc,sc) in scenario_data["Generation"]["Scenarios"];for (k_yr,yr) in sc; for (k_cunt,cuntree) in yr;
        filter!(:Generation_Type=>x->x!="Battery", cuntree);
        filter!(:Generation_Type=>x->x!="VOLL", cuntree);
        push!(cuntree,["SLACK",1000000])
    end;end;end
        push!(scenario_data["Generation"]["keys"],"SLACK")
        delete!(scenario_data["Generation"]["costs"],"VOLL");
        push!(scenario_data["Generation"]["costs"],"SLACK"=>maximum(values(scenario_data["Generation"]["costs"])))

    ####################### Freeze offshore expansion of data #################
    scenario_data=freeze_offshore_expansion(s["nodes"], scenario_data)
    return scenario_data
end

#=
#seperates wfs from genz and defines markets/wfs zones
function main_ACDC_wgen_types(rt_ex,argz, s)
    data, ics_ac, ics_dc, nodes = get_topology_data(rt_ex, s)#topology.m file
    scenario_data = get_scenario_data(nodes)#scenario time series
	###########################################################################
	markets,all_gens,map_gen_types = gen_types(argz["owpp_mva"],nodes,data, scenario_data,s)
    #################### Calculates cable options for AC lines
    data=AC_cable_options(data,argz["candidate_ics_ac"],ics_ac,data["baseMVA"])
    #################### Calculates cable options for DC lines
    data=DC_cable_options(data,argz["candidate_ics_dc"],ics_dc,data["baseMVA"])
    #if (haskey(s, "wf_circuit") && length(s["wf_circuit"])>0);data=keep_only_hm_cables(s,data);end#if home market reduce to only those in 
    additional_params_PMACDC(data)
    print_topology_data_AC(data,markets)#print to verify
    print_topology_data_DC(data,markets)#print to verify
    ##################### load time series data ##############################
    scenario_data = load_time_series_gentypes(rt_ex,argz, scenario_data,markets)
    ##################### multi period setup #################################
	s=update_settings_wgenz(s, argz, data)
    mn_data, xd, map_gen_types  = multi_period_setup_wgen_type(scenario_data, data, all_gens, markets, argz, s, map_gen_types)
	for (g,gen) in mn_data["nw"]
		gen["gen"]=Dict(gen["gen"]);end
	s["xd"]=xd
	push!(s,"max_invest_per_year"=>max_invest_per_year(argz))
    return  mn_data, data, argz, s,map_gen_types
end=#

function freeze_offshore_expansion(nodes, scenario_data)
	developement_zones=filter(:type=>x->x==0,nodes)[!,:country]
	Base_offshore=Dict();for (k_cunt,cuntree) in scenario_data["Generation"]["Scenarios"]["Base"]["2020"]
	push!(Base_offshore,k_cunt=>filter(:Generation_Type=>x->x=="Offshore Wind", cuntree))end;
	for (k_sc,sc) in scenario_data["Generation"]["Scenarios"]
        for (k_yr,yr) in sc;
            for (k_cunt,cuntree) in yr;
                if (issubset([k_cunt],developement_zones))
                    owf=findfirst(x->x=="Offshore Wind", cuntree[!,:Generation_Type]);
                    scenario_data["Generation"]["Scenarios"][k_sc][k_yr][k_cunt][!,:Capacity][owf]=Base_offshore[k_cunt][!,"Capacity"][1]
                end;end;end;end
                return scenario_data
            end

function keep_only_hm_cables(s,data)
    for (b,br) in data["branchdc_ne"]
        if (is_intra_zonal(br["fbusdc"],br["tbusdc"],s["home_market"]))
            #br["rateA"]=br["rateB"]=br["rateC"]=0
            br["status"]=0
        end;end;
    for (b,br) in data["ne_branch"]
        if (is_intra_zonal(br["f_bus"],br["t_bus"],s["home_market"]))
            #br["rateA"]=br["rateB"]=br["rateC"]=0
            br["br_status"]=0
        end;end;
    return data
end

function is_intra_zonal(f,t,zms)
    for zm in zms
        if (issubset([f,t],zm))
            #println(f);println(t)
            return true
        end
    end
    return false
end

function update_settings(s, argz, data)
    s["genz"]=argz["genz"]
    s["wfz"]=argz["wfz"]
    s["ic_lim_onshore"]=argz["conv_lim_onshore"]/data["baseMVA"]
    s["ic_lim_offshore"]=argz["conv_lim_offshore"]/data["baseMVA"]
    s["rad_lim"]=maximum([b["rate_a"] for (k,b) in data["ne_branch"]])
    s["scenarios_length"] = length(argz["scenario_names"])
    s["years_length"] = length(argz["scenario_years"])
    s["hours_length"] = argz["ls"]
    return s
end


function update_settings_wgenz(s, argz, data)
    s["ic_lim_onshore"]=argz["conv_lim_onshore"]/data["baseMVA"]
    s["ic_lim_offshore"]=argz["conv_lim_offshore"]/data["baseMVA"]
    s["rad_lim"]=maximum([b["rate_a"] for (k,b) in data["ne_branch"]])
    s["scenarios_length"] = length(argz["scenario_names"])*length(argz["res_years"])
    s["years_length"] = length(argz["scenario_years"])
    s["hours_length"] = argz["ls"]
    return s
end

#
function update_settings_wgenz(s, data)
    s["ic_lim_onshore"]=s["conv_lim_onshore"]/data["baseMVA"]
    s["ic_lim_offshore"]=s["conv_lim_offshore"]/data["baseMVA"]
    s["rad_lim"]=maximum([b["rate_a"] for (k,b) in data["ne_branch"]])
    s["scenarios_length"] = length(s["scenario_names"])*length(s["res_years"])
    s["years_length"] = length(s["scenario_years"])
    return s
end

########################################## Cables ###############################################
###################### HVAC/HVDC
#loads .m result and filters candidates
function filter_mfile_cables(rt_ex)
    nodes = DataFrames.DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "node_generation")...)
	edges = DataFrames.DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "connections_acdc")...)
    file = rt_ex*"topology.m"
	data = PowerModels.parse_file(file)
	data,ics_ac=filter_AClines(data,edges,nodes)
	data,ics_dc=filter_DClines(data,edges,nodes)
    return data, ics_ac, ics_dc, nodes
end

#utm coordinates from gps
function utm_gps2xy(lla,north_south::Bool=true,zone_utm::Int64=31)
    utm_desired = Geodesy.UTMfromLLA(zone_utm, north_south, Geodesy.wgs84)#sets UTM zone
    utm = utm_desired(Geodesy.LLA(first(lla),last(lla)))#coverts to cartesian
    return utm
end

################################# HVAC ##############################
function filter_AClines(data,edges,nodes)
    #size and length
    ics_ac=Tuple{Int64,Int64}[]
    acc=filter(x->!ismissing(x),edges[!,"AC_mva"])
    for (k, s) in enumerate(acc)
        from_xy=utm_gps2xy((nodes[!,"lat"][edges[!,"AC_from"][k]],nodes[!,"long"][edges[!,"AC_from"][k]]))
        to_xy=utm_gps2xy((nodes[!,"lat"][edges[!,"AC_to"][k]],nodes[!,"long"][edges[!,"AC_to"][k]]))
        push!(ics_ac,(s,round(Int64,Geodesy.euclidean_distance(from_xy, to_xy, 31, true, Geodesy.wgs84)/1000*1.25)))
    end

    #filter ac connections
    accbles2keep_ne=Dict[];accbles2keep=Dict[]
    acc=filter!(x->!ismissing(x),edges[!,"AC_mva"])
    for (r,s) in enumerate(acc)
        for (k,b) in data["ne_branch"]
            if (edges[!,"AC_from"][r]==b["f_bus"] && edges[!,"AC_to"][r]==b["t_bus"])
                push!(accbles2keep_ne,deepcopy(b))
                break;
            end
        end
    end
    data["ne_branch"]=Dict{String,Any}()
    for (k,c) in enumerate(accbles2keep_ne)
        c["source_id"][2]=k
        push!(data["ne_branch"],string(k)=>c)
    end

	for (r,s) in enumerate(acc)
        for (k,b) in data["branch"]
            if (edges[!,"AC_from"][r]==b["f_bus"] && edges[!,"AC_to"][r]==b["t_bus"])
                push!(accbles2keep,deepcopy(b))
                break;
            end
        end
    end
    data["branch"]=Dict{String,Any}()
    for (k,c) in enumerate(accbles2keep)
        c["source_id"][2]=k
        push!(data["branch"],string(k)=>c)
    end
    return data, ics_ac
end

function AC_cable_options(data,candidate_ics_ac,ics_ac,pu)
    z_base_ac=(data["bus"]["1"]["base_kv"])^2/pu
    data=additional_candidatesICS_AC(data,candidate_ics_ac,ics_ac)#adds additional candidates
    for (i,bac) in data["ne_branch"]
    data["ne_branch"][i]=candidateIC_cost_impedance_AC(bac,z_base_ac,pu);end
    data["ne_branch"]=unique_candidateIC_AC(data["ne_branch"])#keep only unique candidates
	temp_cables_ne=Dict{String,Any}()
	temp_cables=Dict{String,Any}()
	cable_pu_costs=Dict{String,Any}()
	cable_pu_r=Dict{String,Any}()
	cable_pu_x=Dict{String,Any}()
    cable_rateA=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(data["ne_branch"]), by=x->parse(Int64,x)))
		last(acb)["source_id"][2]=i
		#last(acb)["br_status"]=0
		push!(temp_cables_ne,string(i)=>last(acb))
		for (j,acb_con) in enumerate(sort(OrderedCollections.OrderedDict(data["branch"]), by=x->parse(Int64,x)))
			if (last(acb)["f_bus"]==last(acb_con)["f_bus"] && last(acb)["t_bus"]==last(acb_con)["t_bus"])
				last(acb_con)["source_id"][2]=j
				push!(temp_cables,string(j)=>last(acb_con))
				if (haskey(cable_pu_costs,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])))
					push!(cable_pu_costs[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["construction_cost"]/(last(acb)["rate_a"]))
					push!(cable_pu_r[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["br_r"])
					push!(cable_pu_x[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["br_x"])
                    push!(cable_rateA[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["rate_a"])
				else
					push!(cable_pu_costs,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["construction_cost"]/(last(acb)["rate_a"])])
					push!(cable_pu_r,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["br_r"]])
					push!(cable_pu_x,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["br_x"]])
                    push!(cable_rateA,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["rate_a"]])
				end
				break
			end
		end
	end
	println(cable_pu_costs)
	data["ne_branch"]=temp_cables_ne
	temp_cables2=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(temp_cables), by=x->parse(Int64,x)))
		trms=length(cable_pu_costs[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		trms_total=sum(b for b in cable_pu_costs[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])

		trms_R=length(cable_pu_r[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		trms_total_R=sum(b for b in cable_pu_r[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])

		trms_X=length(cable_pu_x[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		trms_total_X=sum(b for b in cable_pu_x[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])

		last(acb)["br_r"]=(trms_total_R/trms_R)
		last(acb)["br_x"]=(trms_total_X/trms_X)
		last(acb)["cost"]=(trms_total/trms)
        last(acb)["rateA"]=last(acb)["rateB"]=last(acb)["rateC"]=maximum(b for b in cable_rateA[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		last(acb)["source_id"][2]=i
		push!(temp_cables2,string(i)=>last(acb))
	end
	data["branch"]=temp_cables2
    return data
end

#for each AC candidate capacity an appropriate cable is selected and characteristics stored
function candidateIC_cost_impedance_AC(bac,z_base,s_base)
    cb=_ECO.AC_cbl(bac["rate_a"], bac["length"])
    bac["construction_cost"]=cb.costs.ttl
    bac["br_r"]=((cb.elec.ohm/cb.num)*cb.length)/z_base
    bac["br_x"]=((cb.elec.xl/cb.num)*cb.length)/z_base
    bac["rate_c"]=bac["rate_b"]=bac["rate_a"]=(cb.num*cb.elec.mva)/s_base
    return bac
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

function print_topology_data_AC(data_mip,markets_wfs)
    println("%%%%%%%%%%%%%%%%%%%%%%%%% Nodes %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,n) in enumerate(markets_wfs[1])
        println(string(i)*" - "*n)
    end
    for (i,n) in enumerate(markets_wfs[2])
        println(string(length(markets_wfs[1])+i)*" - OWPP_"*n)
    end
    println("%%%%%%%%%%%%%%%%%%%%%%%% Cables %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,br) in sort(OrderedCollections.OrderedDict(data_mip["ne_branch"]), by=x->parse(Int64,x))
        println(string(i)*": "*string(br["f_bus"])*" - "*string(br["t_bus"])*" MVA: "*string(br["rate_a"])*" Length: "*string(br["length"])*" Cost: "*string(br["construction_cost"])*" Status: "*string(br["br_status"]))
    end
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
end

################################# HVDC Cables
function DC_cable_options(data,candidate_ics_dc,ics_dc,pu)
    z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/pu
    data=additional_candidatesICS_DC(data,candidate_ics_dc,ics_dc)#adds additional candidates
    for (i,bdc) in data["branchdc_ne"]
    data["branchdc_ne"][i]=candidateIC_cost_impedance_DC(bdc,z_base_dc);end
    data["branchdc_ne"]=unique_candidateIC_DC(data["branchdc_ne"])#keep only unique candidates
	temp_cables_ne=Dict{String,Any}()
	temp_cables=Dict{String,Any}()
	cable_pu_costs=Dict{String,Any}()
	cable_pu_r=Dict{String,Any}()
    cable_rateA=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(data["branchdc_ne"]), by=x->parse(Int64,x)))
		last(acb)["source_id"][2]=i
		#last(acb)["status"]=0
		push!(temp_cables_ne,string(i)=>last(acb))
		for (j,acb_con) in enumerate(sort(OrderedCollections.OrderedDict(data["branchdc"]), by=x->parse(Int64,x)))
			if (last(acb)["fbusdc"]==last(acb_con)["fbusdc"] && last(acb)["tbusdc"]==last(acb_con)["tbusdc"])
				last(acb_con)["source_id"][2]=j
				push!(temp_cables,string(j)=>last(acb_con))
				break
			end
		end
		if (haskey(cable_pu_costs,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])))
			push!(cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])],last(acb)["cost"]/(last(acb)["rateA"]/pu))
			push!(cable_pu_r[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])],last(acb)["r"])
            push!(cable_rateA[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])],last(acb)["rateA"])
		else
			push!(cable_pu_costs,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["cost"]/(last(acb)["rateA"]/pu)])
			push!(cable_pu_r,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["r"]])
            push!(cable_rateA,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["rateA"]])
		end
	end
	data["branchdc_ne"]=temp_cables_ne
	temp_cables2=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(temp_cables), by=x->parse(Int64,x)))
		trms=length(cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		trms_total=sum(b for b in cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		#trms_max=minimum(b for b in cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		trms_R=length(cable_pu_r[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		trms_total_R=sum(r for r in cable_pu_r[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		last(acb)["cost"]=trms_total/trms
		#last(acb)["cost"]=trms_max
		last(acb)["r"]=trms_total_R/trms_R
		last(acb)["source_id"][2]=i
        last(acb)["rateA"]=maximum(r for r in cable_rateA[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		push!(temp_cables2,string(i)=>last(acb))
	end
	data["branchdc"]=temp_cables2
    return data
end

function filter_DClines(data,edges,nodes)
    #size and length
    ics_dc=Tuple{Int64,Int64}[]
    dcc=filter(x->!ismissing(x),edges[!,"DC_mva"])
    for (k, s) in enumerate(dcc)
        from_xy=utm_gps2xy((nodes[!,"lat"][edges[!,"DC_from"][k]],nodes[!,"long"][edges[!,"DC_from"][k]]))
        to_xy=utm_gps2xy((nodes[!,"lat"][edges[!,"DC_to"][k]],nodes[!,"long"][edges[!,"DC_to"][k]]))
        push!(ics_dc,(s,round(Int64,Geodesy.euclidean_distance(from_xy, to_xy, 31, true, Geodesy.wgs84)/1000*1.25)))
    end

    #filter dc connections
    dccbles2keep_ne=Dict[]
	dccbles2keep=Dict[]
    for (r,s) in enumerate(dcc)
        for (k,b) in data["branchdc_ne"]
            if (edges[!,"DC_from"][r]==b["fbusdc"] && edges[!,"DC_to"][r]==b["tbusdc"])
                push!(dccbles2keep_ne,deepcopy(b))
                break;
            end
        end
		for (k,b) in data["branchdc"]
            if (edges[!,"DC_from"][r]==b["fbusdc"] && edges[!,"DC_to"][r]==b["tbusdc"])
                push!(dccbles2keep,deepcopy(b))
                break;
            end
        end
    end
    data["branchdc_ne"]=Dict{String,Any}()
    for (k,c) in enumerate(dccbles2keep_ne)
        c["source_id"][2]=k
        push!(data["branchdc_ne"],string(k)=>c)
    end
	data["branchdc"]=Dict{String,Any}()
    for (k,c) in enumerate(dccbles2keep)
        c["source_id"][2]=k
        push!(data["branchdc"],string(k)=>c)
    end
    return data, ics_dc
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

#for each DC candidate capacity an appropriate cable is selected and characteristics stored
function candidateIC_cost_impedance_DC(bdc,z_base)
    cb=_ECO.DC_cbl(bdc["rateA"], bdc["length"])
    #bdc["cost"]=cb.costs.cpx_i+cb.costs.cpx_p
    bdc["cost"]=cb.costs.ttl
    #bdc["r"]=((cb.elec.ohm*10^3/cb.num)*cb.length)/z_base
    bdc["r"]=((cb.elec.ohm/cb.num)*cb.length)/z_base
    bdc["rateC"]=bdc["rateB"]=bdc["rateA"]=cb.num*cb.elec.mva
    return bdc
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


function print_topology_data_DC(data_mip,markets_wfs)
    println("%%%%%%%%%%%%%%%%%%%%%%%%% Nodes %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,n) in enumerate(markets_wfs[1])
        println(string(i)*" - "*n)
    end
    for (i,n) in enumerate(markets_wfs[2])
        println(string(length(markets_wfs[1])+i)*" - OWPP_"*n)
    end
    println("%%%%%%%%%%%%%%%%%% Converters %%%%%%%%%%%%%%%%%%%%")
    for (i,cv) in sort(OrderedCollections.OrderedDict(data_mip["convdc"]), by=x->parse(Int64,x))
            println(string(i)*": "*string(cv["cost"]))
    end;
    println("%%%%%%%%%%%%%%%%%%%%%%%% Cables %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,br) in sort(OrderedCollections.OrderedDict(data_mip["branchdc_ne"]), by=x->parse(Int64,x))
        println(string(i)*": "*string(br["fbusdc"])*" - "*string(br["tbusdc"])*" MVA: "*string(br["rateA"])*" Length: "*string(br["length"])*" cost: "*string(br["cost"])*" status: "*string(br["status"]))
    end
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
end


########################################## Generators ##############################################################

#seperates wfs from genz and defines markets/wfs zones
function genz_n_wfs(owpp_mva,nodes,pu,s)
	s["onshore_nodes"]=[];s["offshore_nodes"]=[];
	infinite_grid=sum(owpp_mva)*3
	s,markets_wfs=divide_onshore_offshore(nodes,s)
    genz=[];wfz=[]
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid/pu));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid/pu));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]/pu));end
    return infinite_grid, genz, wfz, markets_wfs
end

function divide_onshore_offshore(nodes,s)
	s["onshore_nodes"]=[];s["offshore_nodes"]=[];
    markets_wfs=[String[],String[]]#UK,DE,DK must be in same order as .m file gens
    for (k,cuntree) in enumerate(nodes[!,"country"])
        if (nodes[!,"type"][k]>0)
        push!(markets_wfs[1],cuntree);push!(s["onshore_nodes"],k);else
        push!(markets_wfs[2],cuntree);push!(s["offshore_nodes"],k);end
    end
	return s,markets_wfs
end

#
function divide_onshore_offshore(s)
	s["onshore_nodes"]=[];s["offshore_nodes"]=[];
    markets_wfs=[String[],String[]]#UK,DE,DK must be in same order as .m file gens
    for (k,cuntree) in enumerate(s["nodes"][!,"country"])
        if (s["nodes"][!,"type"][k]>0)
        push!(markets_wfs[1],cuntree);push!(s["onshore_nodes"],k);else
        push!(markets_wfs[2],cuntree);push!(s["offshore_nodes"],k);end
    end
	return s, markets_wfs
end

function gen_types(owpp_mva,nodes,data,scenario_data,s)
	s,markets_wfs=divide_onshore_offshore(nodes,s)
	map_gen_types=Dict{String,Any}();push!(map_gen_types,"type"=>Tuple[]);
	base_gens=deepcopy(data["gen"])
	all_gens=Dict{String, Any}()
	push!(all_gens,"onshore"=>Dict())
	for (gen,country) in enumerate(markets_wfs[1])
		push!(all_gens["onshore"],country=>Dict())
		if !(haskey(map_gen_types,"countries"));push!(map_gen_types,"countries"=>Dict());end
		if !(haskey(map_gen_types["countries"],country));push!(map_gen_types["countries"],country=>[]);end
		for (t,type) in enumerate(scenario_data["Generation"]["keys"])
			generator_number=string((gen-1)*length(scenario_data["Generation"]["keys"])+t)
			push!(all_gens["onshore"][country],type=>Dict())
			push!(all_gens["onshore"][country][type],generator_number=>copy(base_gens[string(gen)]))
			all_gens["onshore"][country][type][generator_number]["source_id"]=["gen",parse(Int64,generator_number)]
			all_gens["onshore"][country][type][generator_number]["index"]=parse(Int64,generator_number)
			all_gens["onshore"][country][type][generator_number]["gen_status"]=1
			push!(map_gen_types["type"],(generator_number,type))
			push!(map_gen_types["countries"][country],generator_number)
		end
	end

	push!(all_gens,"offshore"=>Dict())
	for (gen,country) in enumerate(markets_wfs[2])
		if !(haskey(all_gens["offshore"],country))
			push!(all_gens["offshore"],country=>Dict());end
		if !(haskey(map_gen_types,"offshore"));push!(map_gen_types,"offshore"=>Dict());end
		if !(haskey(map_gen_types["offshore"],country));push!(map_gen_types["offshore"],country=>[]);end
		generator_number=string(length(all_gens["onshore"])*length(scenario_data["Generation"]["keys"])+gen)
		push!(all_gens["offshore"][country],generator_number=>copy(base_gens[string(length(markets_wfs[1])+gen)]))
		all_gens["offshore"][country][generator_number]["source_id"]=["gen",parse(Int64,generator_number)]
		all_gens["offshore"][country][generator_number]["index"]=parse(Int64,generator_number)
		all_gens["offshore"][country][generator_number]["gen_status"]=1
		push!(map_gen_types["type"],(generator_number,"Offshore Wind"))
		push!(map_gen_types["offshore"][country],generator_number)
	end
    map_gen_types["markets"]=markets_wfs
	map_gen_types["costs"]=scenario_data["Generation"]["costs"]
    return markets_wfs,all_gens, map_gen_types
end
#
function gen_types(data,scenario_data,s)
	s, markets_wfs = divide_onshore_offshore(s)
	s["map_gen_types"]=Dict{String,Any}();push!(s["map_gen_types"],"type"=>Tuple[]);
	base_gens=deepcopy(data["gen"])
	all_gens=Dict{String, Any}()
	push!(all_gens,"onshore"=>Dict())
	for (gen,country) in enumerate(markets_wfs[1])
		push!(all_gens["onshore"],country=>Dict())
		if !(haskey(s["map_gen_types"],"countries"));push!(s["map_gen_types"],"countries"=>Dict());end
		if !(haskey(s["map_gen_types"]["countries"],country));push!(s["map_gen_types"]["countries"],country=>[]);end
		for (t,type) in enumerate(scenario_data["Generation"]["keys"])
			generator_number=string((gen-1)*length(scenario_data["Generation"]["keys"])+t)
			push!(all_gens["onshore"][country],type=>Dict())
			push!(all_gens["onshore"][country][type],generator_number=>copy(base_gens[string(gen)]))
			all_gens["onshore"][country][type][generator_number]["source_id"]=["gen",parse(Int64,generator_number)]
			all_gens["onshore"][country][type][generator_number]["index"]=parse(Int64,generator_number)
			all_gens["onshore"][country][type][generator_number]["gen_status"]=1
			push!(s["map_gen_types"]["type"],(generator_number,type))
			push!(s["map_gen_types"]["countries"][country],generator_number)
		end
	end

	push!(all_gens,"offshore"=>Dict())
	for (gen,country) in enumerate(markets_wfs[2])
		if !(haskey(all_gens["offshore"],country))
			push!(all_gens["offshore"],country=>Dict());end
		if !(haskey(s["map_gen_types"],"offshore"));push!(s["map_gen_types"],"offshore"=>Dict());end
		if !(haskey(s["map_gen_types"]["offshore"],country));push!(s["map_gen_types"]["offshore"],country=>[]);end
		generator_number=string(length(all_gens["onshore"])*length(scenario_data["Generation"]["keys"])+gen)
		push!(all_gens["offshore"][country],generator_number=>copy(base_gens[string(length(markets_wfs[1])+gen)]))
		all_gens["offshore"][country][generator_number]["source_id"]=["gen",parse(Int64,generator_number)]
		all_gens["offshore"][country][generator_number]["index"]=parse(Int64,generator_number)
		all_gens["offshore"][country][generator_number]["gen_status"]=1
		push!(s["map_gen_types"]["type"],(generator_number,"Offshore Wind"))
		push!(s["map_gen_types"]["offshore"][country],generator_number)
	end
    s["map_gen_types"]["markets"]=markets_wfs
	s["map_gen_types"]["costs"]=scenario_data["Generation"]["costs"]
    return all_gens, s
end

################################################ Converters #####################################
#adds DC grid to PMACDC
function additional_params_PMACDC(data)
    _PMACDC.process_additional_data!(data)#add extra DC model data
    converter_parameters_rxb(data)#sets converter parameters for loss calc
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


############################################# Convex solution to MIP ###############################################
function convex2mip(result_mip, data, mn_data, s)
    s["agent"]=""
    s["relax_problem"]=false
    s["output"]["duals"]=false
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
    return ac_cables
end

function remove_integers(result_mip,mn_data,data,s)
    for (sc,tss) in sort(OrderedCollections.OrderedDict(mn_data["scenario"]), by=x->parse(Int64,x))
        for (t,ts) in sort(OrderedCollections.OrderedDict(tss), by=x->parse(Int64,x))
            #dc cables
            for (bc,brc) in data["branchdc"]
                for (b,br) in result_mip["solution"]["nw"][string(ts)]["branchdc_ne"];
                    if (brc["fbusdc"]==data["branchdc_ne"][b]["fbusdc"] && brc["tbusdc"]==data["branchdc_ne"][b]["tbusdc"])
                        if (br["isbuilt"]>0)
                            s["xd"]["branchdc"][bc]["rateA"][1,ts]=data["branchdc_ne"][b]["rateA"]
                            s["xd"]["branchdc"][bc]["r"][1,ts]=data["branchdc_ne"][b]["r"]
                            s["xd"]["branchdc"][bc]["cost"][1,ts]=0.0;
                            break
                        else
                            s["xd"]["branchdc"][bc]["cost"][1,ts]=0.0;
                            s["xd"]["branchdc"][bc]["rateA"][1,ts]=0.0
                        end
                    end;end;end

            #ac cables
            for (bc,brc) in data["branch"] 
                if (haskey(result_mip["solution"]["nw"][string(ts)],"ne_branch"))
                for (b,br) in result_mip["solution"]["nw"][string(ts)]["ne_branch"];
                        if (brc["f_bus"]==data["ne_branch"][b]["f_bus"] && brc["t_bus"]==data["ne_branch"][b]["t_bus"])                      
                                if (br["built"]>0)
                                    s["xd"]["branch"][bc]["rateA"][1,ts]=data["ne_branch"][b]["rate_a"]
                                    s["xd"]["branch"][bc]["br_r"][1,ts]=data["ne_branch"][b]["br_r"]
                                    s["xd"]["branch"][bc]["cost"][1,ts]=0.0;
                                    break
                                else
                                    s["xd"]["branch"][bc]["rateA"][1,ts]=0
                                    s["xd"]["branch"][bc]["cost"][1,ts]=0.0;
                                end
                        
                    end;end;end;end
    end;end
    return s, mn_data
end



function set_intra_zonal_grid(result_mip,mn_data,s)
    for (sc,tss) in sort(OrderedCollections.OrderedDict(mn_data["scenario"]), by=x->parse(Int64,x))
        for (t,ts) in sort(OrderedCollections.OrderedDict(tss), by=x->parse(Int64,x))
            #dc cables
            for (b,br) in result_mip["solution"]["nw"][string(ts)]["branchdc_ne"];
                if (is_intra_zonal(mn_data["nw"][string(ts)]["branchdc_ne"][b]["fbusdc"],mn_data["nw"][string(ts)]["branchdc_ne"][b]["tbusdc"],s["home_market"]))
                    if (br["isbuilt"]>0)
                            s["xd"]["branchdc_ne"][b]["cost"][1,ts]=0.0;
                    else
                            s["xd"]["branchdc_ne"][b]["cost"][1,ts]=s["xd"]["branchdc_ne"][b]["cost"][1,ts]*100;
                    end
                end;end;

            #ac cables
            if (haskey(result_mip["solution"]["nw"][string(ts)],"ne_branch"))
            for (b,br) in result_mip["solution"]["nw"][string(ts)]["ne_branch"];
                if (is_intra_zonal(mn_data["nw"][string(ts)]["ne_branch"][b]["f_bus"],mn_data["nw"][string(ts)]["ne_branch"][b]["t_bus"],s["home_market"]))
                if (br["built"]>0)
                        s["xd"]["ne_branch"][b]["cost"][1,ts]=0.0;
                else
                    s["xd"]["ne_branch"][b]["construction_cost"][1,ts]=s["xd"]["ne_branch"][b]["construction_cost"][1,ts]*100;
                end;end;
            end;end
            #converters
            if (haskey(result_mip["solution"]["nw"][string(ts)],"convdc"))
                for (c,cnv) in result_mip["solution"]["nw"][string(ts)]["convdc"];
                    if (cnv["p_pacmax"]>0)
                            s["xd"]["convdc"][c]["Pacmin"][1,ts]=cnv["p_pacmax"];
                    end;end;end
    end;end
    return mn_data, s
end

function set_inter_zonal_grid(result_mip,mn_data,s)
    for (sc,tss) in sort(OrderedCollections.OrderedDict(mn_data["scenario"]), by=x->parse(Int64,x))
        for (t,ts) in sort(OrderedCollections.OrderedDict(tss), by=x->parse(Int64,x))
            #dc cables
            for (b,br) in result_mip["solution"]["nw"][string(ts)]["branchdc_ne"];
                if !(is_intra_zonal(mn_data["nw"][string(ts)]["branchdc_ne"][b]["fbusdc"],mn_data["nw"][string(ts)]["branchdc_ne"][b]["tbusdc"],s["home_market"]))
                    if (br["isbuilt"]>0)
                            s["xd"]["branchdc_ne"][b]["cost"][1,ts]=0.0;
                    else
                            s["xd"]["branchdc_ne"][b]["cost"][1,ts]=s["xd"]["branchdc_ne"][b]["cost"][1,ts]*100;
                    end
                end;end;

            #ac cables
            if (haskey(result_mip["solution"]["nw"][string(ts)],"ne_branch"))
            for (b,br) in result_mip["solution"]["nw"][string(ts)]["ne_branch"];
                if !(is_intra_zonal(mn_data["nw"][string(ts)]["ne_branch"][b]["f_bus"],mn_data["nw"][string(ts)]["ne_branch"][b]["t_bus"],s["home_market"]))
                if (br["built"]>0)
                        s["xd"]["ne_branch"][b]["construction_cost"][1,ts]=0.0;
                else
                    s["xd"]["ne_branch"][b]["construction_cost"][1,ts]=s["xd"]["ne_branch"][b]["construction_cost"][1,ts]*100;
                end;end;
            end;end
            #converters
           #= if (haskey(result_mip["solution"]["nw"][string(ts)],"convdc"))
                for (c,cnv) in result_mip["solution"]["nw"][string(ts)]["convdc"];
                    if (cnv["p_pacmax"]>0)
                            s["xd"]["convdc"][c]["Pacmin"][1,ts]=cnv["p_pacmax"];
                    end;end;end=#
    end;end
    return mn_data, s
end


function set_rebalancing_grid(result_mip,mn_data,s)
    for (sc,tss) in sort(OrderedCollections.OrderedDict(mn_data["scenario"]), by=x->parse(Int64,x))
        for (t,ts) in sort(OrderedCollections.OrderedDict(tss), by=x->parse(Int64,x))
            #dc cables
            for (b,br) in result_mip["solution"]["nw"][string(ts)]["branchdc_ne"];
                if (br["isbuilt"]>0)
                        s["xd"]["branchdc_ne"][b]["cost"][1,ts]=0.0;
                    end
                end;

            #ac cables
            if (haskey(result_mip["solution"]["nw"][string(ts)],"ne_branch"))
            for (b,br) in result_mip["solution"]["nw"][string(ts)]["ne_branch"];
                if (br["built"]>0)
                        s["xd"]["ne_branch"][b]["construction_cost"][1,ts]=0.0;
                end;end;end
            #converters
            if (haskey(result_mip["solution"]["nw"][string(ts)],"convdc"))
                for (c,cnv) in result_mip["solution"]["nw"][string(ts)]["convdc"];
                    if (cnv["p_pacmax"]>0)
                            s["xd"]["convdc"][c]["Pacmin"][1,ts]=cnv["p_pacmax"];
                            s["xd"]["convdc"][c]["Pacmax"][1,ts]=cnv["p_pacmax"];
                    else
                            s["xd"]["convdc"][c]["Pacmin"][1,ts]=0;
                            s["xd"]["convdc"][c]["Pacmax"][1,ts]=0;
                    end;end;end
            #storage
            if (haskey(result_mip["solution"]["nw"][string(ts)],"storage"))
                for (b,strg) in result_mip["solution"]["nw"][string(ts)]["storage"];
                    if (strg["e_absmax"]>0)
                            s["xd"]["storage"][b]["pmin"][1,ts]=strg["e_absmax"];
                            s["xd"]["storage"][b]["pmax"][1,ts]=strg["e_absmax"];
                    else
                            s["xd"]["storage"][b]["pmin"][1,ts]=0;
                            s["xd"]["storage"][b]["pmax"][1,ts]=0;
                    end;end;end
            for wf in s["wfz"]
                s["xd"]["gen"][string(first(wf))]["wf_pmax"][1,ts]=result_mip["solution"]["nw"][string(ts)]["gen"][string(first(wf))]["wf_pacmax"];end;
    end;end
    return mn_data, s
end


function combine_solutions(result_mip_hm,result_mip_wf)
    for (n,nw) in result_mip_wf["solution"]["nw"]
        for (b,br) in nw["branchdc_ne"]
            if (br["isbuilt"]==1)
                result_mip_hm["solution"]["nw"][n]["branchdc_ne"][b]=br
            end
        end
        for (b,br) in nw["ne_branch"]
            if (br["built"]==1)
                result_mip_hm["solution"]["nw"][n]["ne_branch"][b]=br
            end
        end
    end
    return result_mip_hm
end   

function hm_market_prices(result_mip, result_mip_hm_prices)
    for (n,nw) in result_mip_hm_prices["solution"]["nw"]
        for (b,bs) in nw["bus"]
            result_mip["solution"]["nw"][n]["bus"][b]["lam_kcl_r"]=deepcopy(bs["lam_kcl_r"])
        end
    end
    return result_mip
end

function set_cable_impedance(data,result_mip)
    last_step=result_mip["solution"]["nw"][string(maximum(parse.(Int64,keys(result_mip["solution"]["nw"]))))]
    #DC cables
    for (b_ne,br_ne) in last_step["branchdc_ne"]
        if (br_ne["isbuilt"]==1)
            for (b,br) in data["branchdc"]
                if (br["fbusdc"]==data["branchdc_ne"][b_ne]["fbusdc"] && br["tbusdc"]==data["branchdc_ne"][b_ne]["tbusdc"])
                    data["branchdc"][b]["r"]=data["branchdc_ne"][b_ne]["r"]
                end
            end
    end;end
    #AC cables
    for (b_ne,br_ne) in last_step["ne_branch"]
        if (br_ne["built"]==1)
            for (b,br) in data["branch"]
                if (br["f_bus"]==data["ne_branch"][b_ne]["f_bus"] && br["t_bus"]==data["ne_branch"][b_ne]["t_bus"])
                    data["branch"][b]["br_r"]=data["ne_branch"][b_ne]["br_r"]
                    data["branch"][b]["br_x"]=data["ne_branch"][b_ne]["br_x"]
                end
            end
    end;end
    return data
end
#seperates wfs from genz and defines markets/wfs zones
function data_update(s,result_mip)
    data, s = get_topology_data(s)#topology.m file
    scenario_data = get_scenario_data(s)#scenario time series
	###########################################################################
	all_gens,s = gen_types(data,scenario_data,s)
    #################### Calculates cable options for AC lines
    data = AC_cable_options(data,s["candidate_ics_ac"],s["ics_ac"],data["baseMVA"])
    #################### Calculates cable options for DC lines
    data = DC_cable_options(data,s["candidate_ics_dc"],s["ics_dc"],data["baseMVA"])
    ############# Sets convex able impedance to the MIP solution ############## 
    data = set_cable_impedance(data, result_mip)
    additional_params_PMACDC(data)
    #print_topology_data_AC(data,s["map_gen_types"]["markets"])#print to verify
    #print_topology_data_DC(data,s["map_gen_types"]["markets"])#print to verify
    ##################### load time series data ##############################
    scenario_data = load_time_series_gentypes(s, scenario_data)
    ##################### multi period setup #################################
	s = update_settings_wgenz(s, data)
    mn_data, s  = multi_period_setup_wgen_type(scenario_data, data, all_gens, s);
	push!(s,"max_invest_per_year"=>max_invest_per_year(s))
    return  mn_data, data, s
end

#seperates wfs from genz and defines markets/wfs zones
function data_setup_nodal(s)
    data, s = get_topology_data(s)#topology.m file
    scenario_data = get_scenario_data(s)#scenario time series
	###########################################################################
	all_gens,s = gen_types(data,scenario_data,s)
    #################### Calculates cable options for AC lines
    data = AC_cable_options(data,s["candidate_ics_ac"],s["ics_ac"],data["baseMVA"])
    #################### Calculates cable options for DC lines
    data = DC_cable_options(data,s["candidate_ics_dc"],s["ics_dc"],data["baseMVA"])
    #if (haskey(s, "wf_circuit") && length(s["wf_circuit"])>0);data=keep_only_hm_cables(s,data);end#if home market reduce to only those in 
    additional_params_PMACDC(data)
    print_topology_data_AC(data,s["map_gen_types"]["markets"])#print to verify
    print_topology_data_DC(data,s["map_gen_types"]["markets"])#print to verify
    ##################### load time series data ##############################
    scenario_data = load_time_series_gentypes(s, scenario_data)
    ##################### multi period setup #################################
	s = update_settings_wgenz(s, data)
    mn_data, s  = multi_period_setup_wgen_type(scenario_data, data, all_gens, s);
	push!(s,"max_invest_per_year"=>max_invest_per_year(s))
    return  mn_data, data, s
end

function data_setup_zonal(s)
    data, s = get_topology_data(s)#topology.m file
    scenario_data = get_scenario_data(s)#scenario time series
	###########################################################################
	all_gens,s = gen_types(data,scenario_data,s)
    #################### Calculates cable options for AC lines
    data = AC_cable_options(data,s["candidate_ics_ac"],s["ics_ac"],data["baseMVA"])
    #################### Calculates cable options for DC lines
    data = DC_cable_options(data,s["candidate_ics_dc"],s["ics_dc"],data["baseMVA"])
    if (haskey(s, "home_market") && length(s["home_market"])>0);data = keep_only_hm_cables(s,data);end#if home market reduce to only those in 
    additional_params_PMACDC(data)
    print_topology_data_AC(data,s["map_gen_types"]["markets"])#print to verify
    print_topology_data_DC(data,s["map_gen_types"]["markets"])#print to verify
    ##################### load time series data ##############################
    scenario_data = load_time_series_gentypes(s, scenario_data)
    ##################### multi period setup #################################
	s = update_settings_wgenz(s, data)
    mn_data, s  = multi_period_setup_wgen_type(scenario_data, data, all_gens, s);
	push!(s,"max_invest_per_year"=>max_invest_per_year(s))
    return  mn_data, data, s
end