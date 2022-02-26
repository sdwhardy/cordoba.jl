    ################## loads external packages ##############################
    using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames
    using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    #import FlexPlan; const _FP = FlexPlan
    #import InfrastructureModels; const _IM = InfrastructureModels
    include("../aux/post_process/functions.jl")

    ################### ENTSO-E scenario selection and rep years/timeline ###############################
    #scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
    scenario_names=["EU17"]

    scenario_years=["2020","2030","2040"]
    #scenario_years=["2020"]

    scenario_planning_horizon=30

    ################### Input data ########################
    owpp_mva=[10000]#mva of wf (in de)
    wf_rad=[(4000,470),(4000,145),(4000,246)]#DE
    #wf_rad=[(4000,160),(4000,45),(4000,602)]#BE

    #interconnectors format: (mva,km)
    ics=[(4000,550),(4000,760),(4000,250)];#UK,DE,DK
    #ics=[(4000,205),(4000,596),(4000,647)];#UK,BE,DK
    ics=vcat(ics,wf_rad);
    conv_lim=4000

    #Countries of nodes
    markets_wfs=[["UK","DE","DK"],["DE"]]#UK,DE,DK must be in same order as .m file gens
    #markets_wfs=[["UK","BE","DK"],["BE"]]#UK,BE,DK must be in same order as .m file gens
    infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

    ##################### load time series data ##############################
    k=5#number of representative days modelled (24 hours per day)
    scenario_data=load("./test/data/input/BE_DE_UK_DK_islands/EUSTDG_17to20_TS_k"*string(k)*"_wBE.jld2")
    #keep only specified scenarios
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],scenario_years));delete!(scenario_data[k],y);end; end;end;end
    #scenario_data["EU17"]["2020"]=scenario_data["EU17"]["2020"][1:108,:]
    #scenario_data["EU19"]["2020"]=scenario_data["EU19"]["2020"][1:2,:]

    ##################### Find minimum length scenario and Make all scenarios the same length
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end


    ################################################# Topology Set up #######################################
    ################## reads .m input file name ######################
    #casename = "test_convex_conv_set"
    casename = "test_convex_conv_set_nippon"
    file = "./test/data/input/$casename.m"
    data_mip = _PM.parse_file(file)#load data in PM format


    z_base_dc=(data_mip["busdc"]["1"]["basekVdc"])^2/data_mip["baseMVA"]
    ######################## seperate gnz and wfs ###################
    genz=[];wfz=[]
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid/data_mip["baseMVA"]));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid/data_mip["baseMVA"]));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]/data_mip["baseMVA"]));end
    #################### Calculates cable options for DC lines
    candidate_ics=[1,4/5,3/5,1/2]#Candidate Cable sizes
    #candidate_ics=[1,4/5]#Candidate Cable sizes
    data_mip=_CBD.additional_candidatesICS_DC(data_mip,candidate_ics,ics)#adds additional candidates
    for (i,bdc) in data_mip["branchdc_ne"]
    data_mip["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance_DC(bdc,z_base_dc);end
    data_mip["branchdc_ne"]=_CBD.unique_candidateIC_DC(data_mip["branchdc_ne"])#keep only unique candidates
    data_mip["branchdc_ne"]=_CBD.cable_reorder_rename(data_mip["branchdc_ne"])#ensure remaining n candidates are ordered 1 to n
    print_topology_data(data_mip,markets_wfs)#print to verify

    ########## remove and rename cables here (case specific adjustments)
    #=delete!(data_mip["branchdc_ne"],"1");delete!(data_mip["branchdc_ne"],"2");delete!(data_mip["branchdc_ne"],"3");
    delete!(data_mip["branchdc_ne"],"4");delete!(data_mip["branchdc_ne"],"5");delete!(data_mip["branchdc_ne"],"6");
    delete!(data_mip["branchdc_ne"],"7");delete!(data_mip["branchdc_ne"],"8");delete!(data_mip["branchdc_ne"],"9");
    delete!(data_mip["branchdc_ne"],"10");delete!(data_mip["branchdc_ne"],"11");delete!(data_mip["branchdc_ne"],"12")
    delete!(data_mip["branchdc_ne"],"13");delete!(data_mip["branchdc_ne"],"14");delete!(data_mip["branchdc_ne"],"15");

    delete!(data_mip["branchdc_ne"],"16");delete!(data_mip["branchdc_ne"],"17");delete!(data_mip["branchdc_ne"],"18");
    delete!(data_mip["branchdc_ne"],"19");delete!(data_mip["branchdc_ne"],"20");delete!(data_mip["branchdc_ne"],"21");
    delete!(data_mip["branchdc_ne"],"22");delete!(data_mip["branchdc_ne"],"23");delete!(data_mip["branchdc_ne"],"24")

    data_mip["branchdc_ne"]["1"]=data_mip["branchdc_ne"]["17"]
    data_mip["branchdc_ne"]["2"]=data_mip["branchdc_ne"]["18"]
    data_mip["branchdc_ne"]["3"]=data_mip["branchdc_ne"]["19"]
    data_mip["branchdc_ne"]["4"]=data_mip["branchdc_ne"]["20"]
    #data_mip["branchdc_ne"]["5"]=data_mip["branchdc_ne"]["17"]
    delete!(data_mip["branchdc_ne"],"17");delete!(data_mip["branchdc_ne"],"18");
    delete!(data_mip["branchdc_ne"],"19");delete!(data_mip["branchdc_ne"],"20");=#
    ##############################################
    #################### Calculates cable options for AC lines
    z_base_ac=(data_mip["bus"]["1"]["base_kv"])^2/data_mip["baseMVA"]
    ac_candidates=[1,4/5,3/5,1/2]#Candidate Cable sizes
    #ac_candidates=[2]#Candidate Cable sizes
    data_mip=_CBD.additional_candidatesICS_AC(data_mip,ac_candidates,wf_rad)#adds additional candidates
    for (i,bac) in data_mip["ne_branch"]
    data_mip["ne_branch"][i]=_CBD.candidateIC_cost_impedance_AC(bac,z_base_ac,data_mip["baseMVA"]);end
    data_mip["ne_branch"]=_CBD.unique_candidateIC_AC(data_mip["ne_branch"])#keep only unique candidates
    data_mip["ne_branch"]=_CBD.cable_reorder_rename(data_mip["ne_branch"])#ensure remaining n candidates are ordered 1 to n
    print_topology_data_AC(data_mip,markets_wfs)#print to verify
    base_price=deepcopy(data_mip)

    _PMACDC.process_additional_data!(data_mip)#add extra DC model data
    _CBD.converter_parameters_rxb(data_mip)#sets converter parameters for loss calc

    #################### Multi-period input parameters #######################
    all_scenario_data,data_mip,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_mip)
    scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
    extradata,data_mip = _CBD.create_profile_sets_mesh(dim, data_mip, all_scenario_data, markets_wfs, infinite_grid, [data_mip["baseMVA"]])

    #################### Scale cost data
    #[println(b*" "*string(br["cost"])) for (b,br) in data_mip["convdc"]];println()
    _CBD.scale_cost_data_2hourly!(data_mip, scenario)#infrastructure investments
    #[println(b*" "*string(br["cost"])) for (b,br) in data_mip["convdc"]];println()
    _CBD.scale_cost_data_2yearly!(extradata, scenario)#energy cost benefits

    # Create data dictionary where time series data is included at the right place
    mn_data_mip = _PMACDC.multinetwork_data(data_mip, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()

    # scale all to NPV
    #mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)#sum of years must equal total
    mn_data_mip= _CBD.npvs_costs_datas_wREZ(mn_data_mip, scenario, scenario_years)#sum of years must equal total
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    mn_data_mip= _CBD.npvs_costs_datas_4mip(mn_data_mip, scenario, scenario_years)#future investment at y scaled to year y=0
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["branchdc_ne"]];println()
    #[println(k*" "*b*" "*string(br["construction_cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["ne_branch"]];println()
    ############################ run optimization

        #################### optimation settings options
        #s = Dict("output" => Dict("branch_flows" => false, "duals"=>false),"fixed_variables" => Dict{String,Any}(),"agent" => "","relax_problem" => false,
        #"conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls,
        #"corridor_limit" => true, "ic_lim"=>conv_lim/100, "strg_lim_offshore"=>0.2, "strg_lim_onshore"=>10, "compare_mode" => false, "dual_update"=>false)

        ################### select solver
        #gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

    s = Dict("output" => Dict("branch_flows" => false, "duals"=>false),"fixed_variables" => Dict{String,Any}(),"agent" => "","relax_problem" => false,
    "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls,
    "corridor_limit" => true, "genz"=>genz, "wfz"=>wfz, "ic_lim"=>conv_lim/100, "strg_lim_offshore"=>0.2, "strg_lim_onshore"=>10, "compare_mode" => false, "dual_update"=>false, "max_invest_per_year"=>10000)
    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
    #
    #result_mip = _CBD.acdc_tnep_convex_conv_strg_npv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#-16164.432804419186 21-13-17
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#-16164.432804419186 21-13-17

    #result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

    #################### fix binaries and converters to solution and solve
    #=s["fixed_variables"] = _CBD.fix_cables_and_converters(mn_data_mip["nw"],s["fixed_variables"], [17], [0,37.41,0,37.41])
    a="fixed_cables_cons"
    s["agent"]=a
    result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function=#
    ####################################################################################
    print_solution_data_DC(result_mip, data_mip)
    print_solution_data_AC(result_mip, data_mip)



#######################################################  testing of formulations below ##################################################################
result_mip["solution"]["nw"][string(1)]
    investment=0
    energy_revenue=0
    for (j,b) in base_price["branchdc_ne"]
        if (result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branchdc_ne"][j]["isbuilt"]>0)
            investment=investment+b["cost"]
        end
    end
    for (j,b) in base_price["ne_branch"]
        if (result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["ne_branch"][j]["built"]>0)
            investment=investment+b["construction_cost"]
        end
    end
    for (j,b) in base_price["convdc"]
            investment=investment+b["cost"]*result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["convdc"][j]["p_pacmax"]
    end
    for (j,b) in base_price["storage"]
            investment=investment+b["cost"]*result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["storage"][j]["e_absmax"]
    end
    investment=investment+base_price["gen"]["4"]["invest"]*result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["gen"]["4"]["wf_pacmax"]

    for (i,nw) in  mn_data_mip["nw"]
        for (j,b) in nw["gen"]
            #println(string(i)*" "*string(j)*" "*string(b["cost"][1])*" "*string(result_mip["solution"]["nw"][i]["gen"][j]["pg"]))
            energy_revenue=energy_revenue+b["cost"][1]*result_mip["solution"]["nw"][i]["gen"][j]["pg"]*mn_data_mip["scenario_prob"]["1"]
        end
    end
    println(string(energy_revenue))
    println(string(investment))
    println(string(energy_revenue+investment))
    println(string(result_mip["objective"]))
