    ################## loads external packages ##############################
    using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames, XLSX
    using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    #########################################################################

    ################### ENTSO-E scenario selection and rep years/timeline ###############################
    scenario_planning_horizon=30
    k=2#number of representative days modelled (24 hours per day)
    #scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
    scenario_names=["EU17"]
    scenario_years=["2020","2030","2040"]
    #scenario_years=["2020"]
    ###########################################################################

    ################# Load topology files ###################################
    rt=pwd()
    include(rt*"\\test\\binaries\\print_m_file.jl")
    include(rt*"\\aux\\post_process\\functions.jl")

    rt_ex=rt*"\\test\\data\\input\\admm_tests\\"
    nodes = DataFrame(XLSX.readtable(rt_ex*"admm_intercons.xlsx", "node_generation")...)
    edges = DataFrame(XLSX.readtable(rt_ex*"admm_intercons.xlsx", "connections_acdc")...)
    ac_cable_df = DataFrame(XLSX.readtable(rt_ex*"admm_intercons.xlsx", "CABLES_AC_SET_UP")...)
    dc_cable_df = DataFrame(XLSX.readtable(rt_ex*"admm_intercons.xlsx", "CABLES_DC_SET_UP")...)
    rem_df = DataFrame(XLSX.readtable(rt_ex*"admm_intercons.xlsx", "REMAINDER")...)
    ppf_mainACDCStorage2mfile(rem_df,ac_cable_df,dc_cable_df,rt_ex)

    casename = "test_case"
    file = rt_ex*"$casename.m"
    data = PowerModels.parse_file(file)
    data,ics_ac=_CBD.filter_AClines(data,edges,nodes)
    data,ics_dc=_CBD.filter_DClines(data,edges,nodes)
    ##########################################################################
    #casename = "test_convex_conv_set_nippon"
    #file = rt_ex*"$casename.m"
    #data = PowerModels.parse_file(file)

    ############### define size and market of genz and wfs ###################
    pu=data["baseMVA"]
    owpp_mva=[4000]#mva of wf (in de)
    conv_lim=4000
    infinite_grid=sum(owpp_mva)*3
    markets_wfs=[String[],String[]]#UK,DE,DK must be in same order as .m file gens
    for (k,cunt) in enumerate(nodes["country"])
        if (nodes["type"][k]>0)
        push!(markets_wfs[1],cunt);else
        push!(markets_wfs[2],cunt);end
    end
    genz=[];wfz=[]
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]));end
    ##########################################################################

    ##################### cables #####################
    #################### Calculates cable options for AC lines
    candidate_ics_ac=[1]
    data=_CBD.additional_candidatesICS_AC(data,candidate_ics_ac,ics_ac)#adds additional candidates
    z_base_ac=(data["bus"]["1"]["base_kv"])^2/pu
    for (i,bac) in data["ne_branch"]
    data["ne_branch"][i]=_CBD.candidateIC_cost_impedance_AC(bac,z_base_ac,pu);end
    data["ne_branch"]=_CBD.unique_candidateIC_AC(data["ne_branch"])#keep only unique candidates
    print_topology_data_AC(data,markets_wfs)#print to verify

    #################### Calculates cable options for DC lines
    candidate_ics_dc=[1]#Candidate Cable sizes
    data=_CBD.additional_candidatesICS_DC(data,candidate_ics_dc,ics_dc)#adds additional candidates
    z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/pu
    for (i,bdc) in data["branchdc_ne"]
    data["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance_DC(bdc,z_base_dc);end
    data["branchdc_ne"]=_CBD.unique_candidateIC_DC(data["branchdc_ne"])#keep only unique candidates
    print_topology_data_DC(data,markets_wfs)#print to verify

    ####################### add DC grid ######################################
    _PMACDC.process_additional_data!(data)#add extra DC model data
    _CBD.converter_parameters_rxb(data)#sets converter parameters for loss calc
    ##########################################################################

    ##################### load time series data ##############################
    scenario_data=load(rt_ex*"EUSTDG_17to20_TS_k"*string(k)*"_wBE.jld2")
    #keep only specified scenarios
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],scenario_years));delete!(scenario_data[k],y);end; end;end;end

    ##################### Find minimum length scenario and Make all scenarios the same length
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end


    #################### Multi-period input parameters #######################
    all_scenario_data,data,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
    extradata,data = _CBD.create_profile_sets_mesh(dim, data, all_scenario_data, markets_wfs, infinite_grid, [data["baseMVA"]])
    #########################################################################

    #################### Scale cost data
    #[println(b*" "*string(br["cost"])) for (b,br) in data_mip["convdc"]];println()
    _CBD.scale_cost_data_2hourly!(data, scenario)#infrastructure investments
    #[println(b*" "*string(br["cost"])) for (b,br) in data_mip["convdc"]];println()
    _CBD.scale_cost_data_2yearly!(extradata, scenario)#energy cost benefits

    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()

    # scale all to NPV
    #mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)#sum of years must equal total
    mn_data= _CBD.npvs_costs_datas_wREZ(mn_data, scenario, scenario_years)#sum of years must equal total
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    mn_data= _CBD.npvs_costs_datas_4mip(mn_data, scenario, scenario_years)#future investment at y scaled to year y=0

    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer)
    s = Dict("output" => Dict("branch_flows" => false, "duals"=>false),"fixed_variables" => Dict{String,Any}(),"agent" => "","relax_problem" => false,
    "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls,
    "corridor_limit" => true, "genz"=>genz, "wfz"=>wfz, "ic_lim"=>conv_lim/pu, "strg_lim_offshore"=>0.2, "strg_lim_onshore"=>10, "compare_mode" => false, "dual_update"=>false, "max_invest_per_year"=>10000)
    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#-16164.432804419186 21-13-17


    print_topology_data_DC(data,[[],[]])
