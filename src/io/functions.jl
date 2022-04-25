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

#seperates wfs from genz and defines markets/wfs zones
function main_ACDC_wgen_types(rt_ex,argz, s)
    ################# Load topology files ###################################
    topology_df(rt_ex, s["relax_problem"], s["AC"])#creates .m file
    data, ics_ac, ics_dc, nodes = filter_mfile_cables(rt_ex)#loads resulting topology and filters for candidate cables
    ############### defines size and market of genz and wfs ###################
	scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4UKBEDEDK.jld2")
	######## Batteries are removed and modeled seperately time series #########
	for (k_sc,sc) in scenario_data["Generation"]["Scenarios"]
		for (k_cunt,cuntree) in sc;filter!(:Generation_Type=>x->x!="Battery", cuntree);end;end
	####################### Freeze offshore expansion of data #################
	developement_zones=filter(:type=>x->x==0,nodes)[!,:country]
	Base_offshore=Dict();for (k_cunt,cuntree) in scenario_data["Generation"]["Scenarios"]["NT2025"]
	push!(Base_offshore,k_cunt=>filter(:Generation_Type=>x->x=="Offshore Wind", cuntree))end;
	for (k_sc,sc) in scenario_data["Generation"]["Scenarios"]
		for (k_cunt,cuntree) in sc;
			if (issubset([k_cunt],developement_zones))
				owf=findfirst(x->x=="Offshore Wind", cuntree[!,:Generation_Type]);
				scenario_data["Generation"]["Scenarios"][k_sc][k_cunt][!,:Capacity][owf]=Base_offshore[k_cunt][!,"Capacity"][1]
		end;end;end
	###########################################################################
	markets,all_gens,map_gen_types = gen_types(argz["owpp_mva"],nodes,data, scenario_data,s)
	map_gen_types["markets"]=markets
	map_gen_types["costs"]=scenario_data["Generation"]["costs"]
    #################### Calculates cable options for AC lines
    data=AC_cable_options(data,argz["candidate_ics_ac"],ics_ac,data["baseMVA"])
    print_topology_data_AC(data,markets)#print to verify
    #################### Calculates cable options for DC lines
    data=DC_cable_options(data,argz["candidate_ics_dc"],ics_dc,data["baseMVA"])
    additional_params_PMACDC(data)
    print_topology_data_DC(data,markets)#print to verify
    ##################### load time series data ##############################
    scenario_data = load_time_series_gentypes(rt_ex,argz, scenario_data,markets)
    ##################### multi period setup #################################
	s["ic_lim"]=argz["conv_lim"]/data["baseMVA"]
    s["rad_lim"]=maximum([b["rate_a"] for (k,b) in data["ne_branch"]])
    s["scenarios_length"] = length(argz["scenario_names"])
    s["years_length"] = length(argz["scenario_years"])
    s["hours_length"] = argz["ls"]
    mn_data, xd, map_gen_types  = multi_period_setup_wgen_type(scenario_data, data, all_gens, markets, argz, s, map_gen_types)
	for (g,gen) in mn_data["nw"]
		gen["gen"]=Dict(gen["gen"]);end
	s["xd"]=xd
	push!(s,"max_invest_per_year"=>max_invest_per_year(argz))
    return  mn_data, data, argz, s,map_gen_types
end

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
				else
					push!(cable_pu_costs,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["construction_cost"]/(last(acb)["rate_a"])])
					push!(cable_pu_r,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["br_r"]])
					push!(cable_pu_x,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["br_x"]])
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
        println(string(i)*": "*string(br["f_bus"])*" - "*string(br["t_bus"])*" MVA: "*string(br["rate_a"])*" Length: "*string(br["length"])*" Cost: "*string(br["construction_cost"]))
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
		else
			push!(cable_pu_costs,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["cost"]/(last(acb)["rateA"]/pu)])
			push!(cable_pu_r,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["r"]])
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
        println(string(i)*": "*string(br["fbusdc"])*" - "*string(br["tbusdc"])*" MVA: "*string(br["rateA"])*" Length: "*string(br["length"])*" cost: "*string(br["cost"]))
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

    return markets_wfs,all_gens, map_gen_types
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

####################################### Print solution ################################

function print_solution_data(result_mip, data_mip, argz)
    println("Description: test-"*string(argz["test"])*" k-"*string(argz["k"])*" years-"*string(argz["scenario_years"])*" scenarios-"*string(argz["scenario_names"]))
    if (haskey(result_mip["solution"]["nw"]["1"],"branch"))
        println("%%%%%%% CONVEX SOLUTION %%%%%%%")
        print_branch(result_mip,argz,data_mip)
        print_branchdc(result_mip,argz,data_mip)
    else
        println("%%%%%%% MIP SOLUTION %%%%%%%")
        print_branch_ne(result_mip,argz,data_mip)
        print_branchdc_ne(result_mip,argz,data_mip)
    end
    print_owpps(result_mip,argz)
    print_converters(result_mip,argz)
    print_storage(result_mip,argz)
    println("objective: "*string(result_mip["objective"])*" achieved in: "*string(result_mip["solve_time"]))
end

function print_branch(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"branch"))
        println("%%%% Cables HVAC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["branch"]), by=x->parse(Int64,x))
            println(string(i)*": "*string(data_mip["branch"][i]["f_bus"])*" - "*string(data_mip["branch"][i]["t_bus"])*" MVA: "*string(br["p_rateAC"]))
        end
        println("%%%% Cables HVAC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["branch"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branch"][i]["f_bus"])*" - "*string(data_mip["branch"][i]["t_bus"])*" MVA: "*string(br["p_rateAC"]))
        end
        println("%%%% Cables HVAC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branch"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branch"][i]["f_bus"])*" - "*string(data_mip["branch"][i]["t_bus"])*" MVA: "*string(br["p_rateAC"]))
        end
    end
end

function print_branchdc(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"branchdc"))
        println("%%%% Cables HVDC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["branchdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branchdc"][i]["fbusdc"])*" - "*string(data_mip["branchdc"][i]["tbusdc"])*" MVA: "*string(br["p_rateA"]))
        end
        println("%%%% Cables HVDC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["branchdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branchdc"][i]["fbusdc"])*" - "*string(data_mip["branchdc"][i]["tbusdc"])*" MVA: "*string(br["p_rateA"]))
        end
        println("%%%% Cables HVDC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branchdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branchdc"][i]["fbusdc"])*" - "*string(data_mip["branchdc"][i]["tbusdc"])*" MVA: "*string(br["p_rateA"]))
        end
    end
end

function print_branch_ne(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"ne_branch"))
        println("%%%% Cables HVAC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["ne_branch"]), by=x->parse(Int64,x))
            if (br["built"]==1)
                println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
        end;end
        println("%%%% Cables HVAC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["ne_branch"]), by=x->parse(Int64,x))
            if (br["built"]==1)
                println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
        end;end
        println("%%%% Cables HVAC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["ne_branch"]), by=x->parse(Int64,x))
            if (br["built"]==1)
                println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
        end;end
    end
end

function print_branchdc_ne(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"branchdc_ne"))
        println("%%%% Cables HVDC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["branchdc_ne"]), by=x->parse(Int64,x))
            if (br["isbuilt"]==1)
                println(string(i)*": "*string(data_mip["branchdc_ne"][i]["fbusdc"])*" - "*string(data_mip["branchdc_ne"][i]["tbusdc"])*" MVA: "*string(data_mip["branchdc_ne"][i]["rateA"])*" cost: "*string(data_mip["branchdc_ne"][i]["cost"]))
        end;end
        println("%%%% Cables HVDC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["branchdc_ne"]), by=x->parse(Int64,x))
            if (br["isbuilt"]==1)
                println(string(i)*": "*string(data_mip["branchdc_ne"][i]["fbusdc"])*" - "*string(data_mip["branchdc_ne"][i]["tbusdc"])*" MVA: "*string(data_mip["branchdc_ne"][i]["rateA"])*" cost: "*string(data_mip["branchdc_ne"][i]["cost"]))
        end;end
        println("%%%% Cables HVDC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branchdc_ne"]), by=x->parse(Int64,x))
            if (br["isbuilt"]==1)
                println(string(i)*": "*string(data_mip["branchdc_ne"][i]["fbusdc"])*" - "*string(data_mip["branchdc_ne"][i]["tbusdc"])*" MVA: "*string(data_mip["branchdc_ne"][i]["rateA"])*" cost: "*string(data_mip["branchdc_ne"][i]["cost"]))
        end;end
    end
end

function print_converters(result_mip,argz)
    if (haskey(result_mip["solution"]["nw"]["1"],"convdc"))
        println("%%%% Converters t0 %%%%")
        for (i,cv) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["convdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(cv["p_pacmax"]))
        end;
        println("%%%% Converters t2 %%%%")
        for (i,cv) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["convdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(cv["p_pacmax"]))
        end;
        println("%%%% Converters tinf %%%%")
        for (i,cv) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["convdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(cv["p_pacmax"]))
        end;
    end
end

function print_storage(result_mip,argz)
    if (haskey(result_mip["solution"]["nw"]["1"],"storage"))
        println("%%%% Storage t0 %%%%")
        for (i,s) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["storage"]), by=x->parse(Int64,x))
                println(string(i)*": "*" MWh: "*string(s["e_absmax"]))
        end
        println("%%%% Storage t2 %%%%")
        for (i,s) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["storage"]), by=x->parse(Int64,x))
                println(string(i)*": "*" MWh: "*string(s["e_absmax"]))
        end
        println("%%%% Storage tinf %%%%")
        for (i,s) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["storage"]), by=x->parse(Int64,x))
                println(string(i)*": "*" MWh: "*string(s["e_absmax"]))
        end
    end
end

function print_owpps(result_mip,argz)
    if (haskey(result_mip["solution"]["nw"]["1"],"gen"))
        println("%%%% OWPPS T0 %%%%")
        for (i,wf) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["gen"]), by=x->parse(Int64,x))
            if (haskey(wf,"wf_pacmax"))
                println(string(i)*": "*string(wf["wf_pacmax"]))
        end;end
        println("%%%% OWPPS T2 %%%%")
        for (i,wf) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["gen"]), by=x->parse(Int64,x))
            if (haskey(wf,"wf_pacmax"))
                println(string(i)*": "*string(wf["wf_pacmax"]))
        end;end
        println("%%%% OWPPS Tinf %%%%")
        for (i,wf) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["gen"]), by=x->parse(Int64,x))
            if (haskey(wf,"wf_pacmax"))
                println(string(i)*": "*string(wf["wf_pacmax"]))
        end;end
    end
end


####################################### Genration/Demand Summary ###################################################

function summarize_generator_solution_data(result_mip, data, argz, map_gen_types,s)#print solution
	gen_tbls=build_generator_tables(result_mip, data)
	gen_by_market=sort_by_country(gen_tbls,map_gen_types,s)
	gen_by_offshore=sort_by_offshore(gen_tbls,map_gen_types,s)
	load_by_market=sort_load_by_country(gen_tbls,map_gen_types)
	gen_consume=Dict()
	push!(gen_consume,"onshore_generation"=>gen_by_market)
	push!(gen_consume,"offshore_generation"=>gen_by_offshore)
	push!(gen_consume,"onshore_demand"=>load_by_market)
	return gen_consume
end
function sort_load_by_country(gen_tbls,map_gen_types)
	per_market=Dict()
	for (s,sc) in gen_tbls
		if !(haskey(per_market,s));push!(per_market,s=>Dict());end
		col_names=names(sc)
		for (cuntree,gens) in map_gen_types["loads"]
			for gen in gens
			if !(haskey(per_market[s],cuntree));push!(per_market[s],cuntree=>DataFrames.DataFrame(Symbol(col_names[1])=>sc[!,Symbol(col_names[1])]));end
				col_num=findfirst(x->x==string(gen),col_names)
				if !(isnothing(col_num))
				per_market[s][cuntree]=hcat(per_market[s][cuntree],DataFrames.DataFrame(Symbol(col_names[col_num])=>sc[!,Symbol(col_names[col_num])]))
			end;end
		end
	end

	return per_market
end

function sort_by_offshore(gen_tbls,map_gen_types,set)
	per_market=Dict()
	for (s,sc) in gen_tbls
		if !(haskey(per_market,s));push!(per_market,s=>Dict());end
		col_names=names(sc)
		for (cuntree,gens) in map_gen_types["offshore"]
			for gen in gens
			if !(haskey(per_market[s],cuntree));push!(per_market[s],cuntree=>DataFrames.DataFrame(Symbol(col_names[1])=>sc[!,Symbol(col_names[1])]));end
				col_num=findfirst(x->x==gen,col_names)
				if !(isnothing(col_num))
				per_market[s][cuntree]=hcat(per_market[s][cuntree],DataFrames.DataFrame(Symbol("Offshore Wind")=>sc[!,Symbol(col_names[col_num])]))
			end;end
		end

		for cuntree_num in set["offshore_nodes"]
			cuntree=map_gen_types["markets"][2][cuntree_num-length(set["onshore_nodes"])]
			if !(haskey(per_market[s],cuntree));push!(per_market[s],cuntree=>DataFrames.DataFrame(Symbol(col_names[1])=>sc[!,Symbol(col_names[1])]));end
			per_market[s][cuntree]=hcat(per_market[s][cuntree],DataFrames.DataFrame(Symbol("Battery")=>sc[!,Symbol("Battery "*string(cuntree_num))]))
		end
	end
	return per_market
end

function sort_by_country(gen_tbls,map_gen_types,set)
	per_market=Dict()
	for (s,sc) in gen_tbls
		if !(haskey(per_market,s));push!(per_market,s=>Dict());end
		col_names=names(sc)
		for (cuntree,gens) in map_gen_types["countries"]
			for gen in gens
			if !(haskey(per_market[s],cuntree));push!(per_market[s],cuntree=>DataFrames.DataFrame(Symbol(col_names[1])=>sc[!,Symbol(col_names[1])]));end
				col_num=findfirst(x->x==gen,col_names)
				if !(isnothing(col_num))
				per_market[s][cuntree]=hcat(per_market[s][cuntree],DataFrames.DataFrame(Symbol(col_names[col_num])=>sc[!,Symbol(col_names[col_num])]))
			end;end
			for (num,type) in map_gen_types["type"];
				cunt_col_names=names(per_market[s][cuntree])
				location=findfirst(x->x==num,cunt_col_names)
				if !(isnothing(location))
					DataFrames.rename!(per_market[s][cuntree], cunt_col_names[location]=>type)
				end
			end
		end

		for cuntree_num in set["onshore_nodes"]
			println(cuntree_num)
			cuntree=map_gen_types["markets"][1][cuntree_num]
			if !(haskey(per_market[s],cuntree));push!(per_market[s],cuntree=>DataFrames.DataFrame(Symbol(col_names[1])=>sc[!,Symbol(col_names[1])]));end
				per_market[s][cuntree]=hcat(per_market[s][cuntree],DataFrames.DataFrame(Symbol("Battery")=>sc[!,Symbol("Battery "*string(cuntree_num))]))
		end
	end
	return per_market
end

function build_generator_tables(result_mip, data)
	gen_per_scenario=Dict{}()
	for (s,sc) in data["scenario"]
		titles=["ts"]
		for (g,gen) in sort!(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["gen"]), by=x->parse(Int64,x));push!(titles,string(g));end
		for (s,bat) in sort!(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["storage"]), by=x->parse(Int64,x));push!(titles,"Battery "*string(s));end
		for (cts,ts) in sort!(OrderedCollections.OrderedDict(sc), by=x->parse(Int64,x))
			ts_row=[];push!(ts_row,ts)
			for (g,gen) in sort!(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(ts)]["gen"]), by=x->parse(Int64,x))
				push!(ts_row,gen["pg"])
			end
			for (s,bat) in sort!(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(ts)]["storage"]), by=x->parse(Int64,x))
				push!(ts_row,-1*bat["ps"])
			end
			titles=hcat(titles,ts_row)
		end

		namelist = Symbol.(titles[1:end,1])
		df  = DataFrames.DataFrame()
		for (i, name) in enumerate(namelist)
    	df[!,name] =  [titles[i,j] for j in 2:length(titles[i,:])];end
		push!(gen_per_scenario,s=>df)
	end
	return gen_per_scenario
end

######################### plotting generation types ###################################

function plot_marginal_price(gen,map_gen_types, country)
    col_names=names(gen[!,2:end])
    for col in col_names;
        if (isapprox(sum(gen[!,col]),0,atol=1))
            DataFrames.select!(gen, DataFrames.Not(Symbol(col)))
    end;end
    marginal_prices=[]
    for row in eachrow(gen)
	    active_gen_types=[]
	    for k in keys(row[2:end])
		    if (row[k]>0)
			    push!(active_gen_types,map_gen_types["costs"][string(k)])
		    end
	    end
	    push!(marginal_prices,maximum(active_gen_types))
    end
    low_rng=minimum(marginal_prices)
    high_rng=maximum(marginal_prices)
    scatter_vec=[
	PlotlyJS.scatter(
	    x=gen[!,:ts], y=marginal_prices,
	    name="Marginal Price", mode="lines",
	    line=PlotlyJS.attr(width=1, color="black")
	)]

	PlotlyJS.plot(
	scatter_vec, PlotlyJS.Layout(yaxis_range=(low_rng, high_rng),yaxis_title="€/MWh",xaxis_title="time steps",title=country))
end


function plot_generation_profile(gen, con, country)
    clrs=generation_color_map()
    col_names=names(gen[!,2:end])
    for col in col_names;
        if (isapprox(sum(gen[!,col]),0,atol=1))
            DataFrames.select!(gen, DataFrames.Not(Symbol(col)))
    end;end
	#battery energy
	con_sum=DataFrames.DataFrame();con_sum[!,:ts]=con[!,:ts]
	if (hasproperty(gen,:Battery))
		bat=gen[!,"Battery"]
		bat_d=[imp>=0 ? imp : 0 for imp in bat]
		if (sum(bat_d)>0);gen[!,"Battery Discharge"]=bat_d;end
		bat_c=[exp<=0 ? exp : 0 for exp in bat]
		if (sum(bat_c)<0);con[!,"Battery Charge"]=bat_c;end
		if (sum(bat_c)<0);con_sum[!,"Battery Charge"]=bat_c;end
		gen=DataFrames.select!(gen,DataFrames.Not(Symbol("Battery")))
	end


    col_names_con=names(con[!,1:end])
    all_con=length(col_names_con)>1 ? abs.(sum(eachcol(con[!,2:end]))) : zeros(Int8,length(con[!,:ts]))
    all_gen=abs.(sum(eachcol(gen[!,2:end])))
    #imported energy
    import_export=all_con.-all_gen
    imp=[imp>=0 ? imp : 0 for imp in import_export]
    if (sum(imp)>0);gen[!,"Import"]=imp;end
    #Set range
    all_gen=abs.(sum(eachcol(gen[!,2:end])))
    rng_gen=maximum(all_gen)
    #exported energy
    exp=[exp<=0 ? exp : 0 for exp in import_export]
    if (sum(exp)<0);con_sum[!,"Export"]=exp;end
	if (sum(all_con)>0);con_sum[!,"Demand"]=-1*all_con;end

	all_con=abs.(sum(eachcol(con_sum[!,2:end])))
    rng_con=maximum(all_con)

    col_names_gen=names(gen[!,2:end])
    col_names_con=names(con_sum[!,2:end])
    scatter_vec_gen=[
        PlotlyJS.scatter(
            x=gen[!,:ts], y=gen[!,Symbol(nm)],
            stackgroup="one", name=String(nm), mode="lines", hoverinfo="x+y",
            line=PlotlyJS.attr(width=0.5, color=clrs[nm])
        ) for nm in col_names_gen]

        scatter_vec_con=[
            PlotlyJS.scatter(
                x=con_sum[!,:ts], y=con_sum[!,Symbol(nm)],
                stackgroup="two", name=String(nm), mode="lines", hoverinfo="x+y",
                line=PlotlyJS.attr(width=0.5, color=clrs[nm])
            ) for nm in col_names_con]

         scatter_vec=vcat(scatter_vec_con,scatter_vec_gen)
            PlotlyJS.plot(
            scatter_vec, PlotlyJS.Layout(yaxis_range=(-1*rng_con, rng_gen),yaxis_title="MW",xaxis_title="time steps",title=country))
end

function generation_color_map()
        color_dict=Dict("Offshore Wind"=>"darkgreen",
        "Onshore Wind"=>"forestgreen",
        "Solar PV"=>"yellow",
        "Solar Thermal"=>"orange",
        "Gas CCGT new"=>"chocolate",
		"Gas OCGT new"=>"chocolate2",
        "Gas CCGT old 1"=>"chocolate3",
        "Gas CCGT old 2"=>"tan1",
        "Gas CCGT present 1"=>"tan2",
        "Gas CCGT present 2"=>"sienna",
        "Reservoir"=>"blue1",
        "Run-of-River"=>"navy",
        "Nuclear"=>"gray69",
        "Other RES"=>"yellowgreen",
        "Gas CCGT new CCS"=>"sienna1",
        "Gas CCGT present 1 CCS"=>"sienna2",
        "Gas CCGT present 2 CCS"=>"sienna3",
        "Battery Discharge"=>"azure",
		"Battery Charge"=>"white",
        "Gas CCGT CCS"=>"sienna4",
        "Demand"=>"grey",
        "Import"=>"red",
        "Export"=>"black")
        return color_dict
end


############################################# Convex solution to MIP ###############################################
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
