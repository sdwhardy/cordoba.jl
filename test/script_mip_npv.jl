    ################## loads external packages ##############################
    using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames
    using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    #import FlexPlan; const _FP = FlexPlan
    #import InfrastructureModels; const _IM = InfrastructureModels
    include("../aux/post_process/functions.jl")

    ################### ENTSO-E scenario description ####################################
    scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
    #scenario_names=["EU17","EU18","EU19","EU20"]
    #scenario_names=["ST18"]

    scenario_years=["2020","2030","2040"]
    #scenario_years=["2040"]
    #scenario_years=["2020"]

    scenario_planning_horizon=30

    ################### Input data for script_test.jl ########################
    owpp_mva=[4000]#mva of wf (in de)
    #owpp_mva=[0]#mva of wf (in de)

    #interconnectors format: (mva,km)
    ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];#UK,DE,DK
    #ics=[(4000,205),(4000,596),(4000,647),(4000,160),(4000,45),(4000,602)];#UK,BE,DK
    #ics=[(4000,145)];
    conv_lim=4000

    #location of nodes
    markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
    #markets_wfs=[["UK","BE","DK"],["BE"]]#must be in same order as .m file gens
    #markets_wfs=[["DE"],["DE"]]
    infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

    ##################### load time series data ##############################
    k=2
    scenario_data=load("./test/data/input/BE_DE_UK_DK_islands/EUSTDG_17to20_TS_k"*string(k)*"_wBE.jld2")
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],scenario_years));delete!(scenario_data[k],y);end; end;end;end
    #d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);end;end

    ##################### Find minimum length scenario
    #scenario_data=load("./test/data/input/EUSTDG_TS_k5.jld2")
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    #scenario_data[_sc][_yr]=scenario_data[_sc][_yr][1:2,:]
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    ##################### Make all scenarios the same length
    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end
    #save("./test/data/input/EUSTDG_TS_k"*string(k)*".jld2",scenario_data)

    #optimation settings
    #s = Dict("output" => Dict("branch_flows" => false), "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls, "corridor_limit" => true, "ic_lim"=>conv_lim/100, "compare_mode" => true)
    s = Dict("output" => Dict("branch_flows" => false, "duals"=>false),"fixed_variables" => Dict{String,Any}(),"agent" => "","relax_problem" => false,
    "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls,
    "corridor_limit" => true, "ic_lim"=>conv_lim/100, "strg_lim_offshore"=>0.2, "strg_lim_onshore"=>10, "compare_mode" => false, "dual_update"=>false)

    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)


    ################################################# MIP TNEP #######################################
    ################## reads .m input file name ######################
    #casename = "test_convex_wf"
    #casename = "test_convset_wf"
    #casename = "test_convex_conv"
    casename = "test_convex_conv_set"
    file = "./test/data/input/$casename.m"
    data_mip = _PM.parse_file(file)#load data in PM format

    #################### Calculates cables for DC lines
    z_base_dc=(data_mip["busdc"]["1"]["basekVdc"])^2/data_mip["baseMVA"]
    #candidate_ics=[1]#Candidate Cable sizes
    #candidate_ics=[4/5]#Candidate Cable sizes
    #candidate_ics=[1,4/5,1/2,1/4]#Candidate Cable sizes
    candidate_ics=[1,4/5,3/5,1/2]#Candidate Cable sizes
    #candidate_ics=[1]#Candidate Cable sizes
    data_mip=_CBD.additional_candidatesICS(data_mip,candidate_ics,ics)#adds additional candidates
    for (i,bdc) in data_mip["branchdc_ne"]
    data_mip["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance(bdc,z_base_dc);end
    data_mip["branchdc_ne"]=_CBD.unique_candidateIC(data_mip["branchdc_ne"])#keep only unique candidates
    print_topology_data(data_mip,markets_wfs)

    ########## remove and rename cables here
    #=delete!(data_mip["branchdc_ne"],"1");delete!(data_mip["branchdc_ne"],"2");
    delete!(data_mip["branchdc_ne"],"6");delete!(data_mip["branchdc_ne"],"3");delete!(data_mip["branchdc_ne"],"4");
    delete!(data_mip["branchdc_ne"],"5");
    delete!(data_mip["branchdc_ne"],"7");delete!(data_mip["branchdc_ne"],"8");
    delete!(data_mip["branchdc_ne"],"9");delete!(data_mip["branchdc_ne"],"13")
    delete!(data_mip["branchdc_ne"],"10");delete!(data_mip["branchdc_ne"],"11");delete!(data_mip["branchdc_ne"],"12")
    delete!(data_mip["branchdc_ne"],"14");delete!(data_mip["branchdc_ne"],"15");delete!(data_mip["branchdc_ne"],"16")
    #delete!(data_mip["branchdc_ne"],"19");delete!(data_mip["branchdc_ne"],"18");delete!(data_mip["branchdc_ne"],"20")
    delete!(data_mip["branchdc_ne"],"21");delete!(data_mip["branchdc_ne"],"23");delete!(data_mip["branchdc_ne"],"22");
    delete!(data_mip["branchdc_ne"],"24")
    data_mip["branchdc_ne"]["1"]=data_mip["branchdc_ne"]["17"]
    data_mip["branchdc_ne"]["2"]=data_mip["branchdc_ne"]["18"]
    data_mip["branchdc_ne"]["3"]=data_mip["branchdc_ne"]["19"]
    data_mip["branchdc_ne"]["4"]=data_mip["branchdc_ne"]["20"]
    #data_mip["branchdc_ne"]["5"]=data_mip["branchdc_ne"]["17"]
    delete!(data_mip["branchdc_ne"],"17");delete!(data_mip["branchdc_ne"],"18");
    delete!(data_mip["branchdc_ne"],"19");delete!(data_mip["branchdc_ne"],"20");=#
    ##############################################

    _PMACDC.process_additional_data!(data_mip)#add extra DC model data
    _CBD.converter_parameters_rxb(data_mip)

    #################### Multi-period input parameters #######################
    all_scenario_data,data_mip,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_mip)
    scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
    extradata,data_mip = _CBD.create_profile_sets_mesh(dim, data_mip, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)

    #################### Scale cost data
    #_CBD.scale_cost_data_2hourly!(data_mip, scenario)
    #_CBD.scale_cost_data_2hourly!(extradata, scenario)
    _CBD.scale_cost_data_2hourly!(data_mip, scenario)
    _CBD.scale_cost_data_2yearlyhourly!(extradata, scenario)

    # Create data dictionary where time series data is included at the right place
    mn_data_mip = _PMACDC.multinetwork_data(data_mip, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))

    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    mn_data_mip= _CBD.npvs_costs_datas_4mip(mn_data_mip, scenario, scenario_years)
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    #run optimization
    #acdc_tnep_convex_conv_strg_npv
    result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

    #################### fix binaries and converters to solution and solve
    s["fixed_variables"] = _CBD.fix_cables_and_converters(mn_data_mip["nw"],s["fixed_variables"], [13,17,21], [40,40,40,40])
    a="fixed_cables_cons"
    s["agent"]=a
    result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
    ###############################################################################
    print_solution_data(result_mip, data_mip)



    ############ testing of formulations
    #result_mip = _CBD.cordoba_mp_acdctnepopf_stoch(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    #result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.LPACCPowerModel, gurobi, multinetwork=true; setting = s)
    #result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.SOCWRPowerModel, gurobi, multinetwork=true; setting = s)

    #print results
    #=println()
    max_key=maximum([parse(Int64,k) for (k,nw) in result_mip["solution"]["nw"]])
    for (k,b) in result_mip["solution"]["nw"][string(max_key)]["branchdc_ne"]
        if (b["isbuilt"]==1)
            println(k*" "*string(mn_data_mip["nw"][string(max_key)]["branchdc_ne"][k]["fbusdc"])*" - "*string(mn_data_mip["nw"][string(max_key)]["branchdc_ne"][k]["tbusdc"])*" MVA: "*string(mn_data_mip["nw"][string(max_key)]["branchdc_ne"][k]["rateA"]*data_mip["baseMVA"]))
            end
    end
    println()
    for (k,b) in result_mip["solution"]["nw"]["1"]["convdc"]
            println(k*" "*" MVA: "*string(b["p_pacmax"]*data_mip["baseMVA"]))=#
    #end
    #[println(string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["fbusdc"])*" - "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["tbusdc"])*" MVA: "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["rateA"]*100)*" Cost: "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["cost"]*scenario["hours"])) for (k,b) in result_mip["solution"]["nw"]["1"]["branchdc_ne"]]

    #nw=mn_data_mip["nw"][i]["branchdc_ne"][j]["tbusdc"]
    revenue=0
    wnd_energy=0
    wnd_revenue_min=0
    wnd_revenue_hm=0
    wnd_revenue_max=0
    capex_cable=0
    capex_conv=0
    for (i,nw) in  mn_data_mip["nw"]
        for (j,b) in nw["branchdc_ne"]
            if (result_mip["solution"]["nw"][i]["branchdc_ne"][j]["isbuilt"]>0)
                #println("nw "*string(i)*" b "*string(j)*" cost "*string(nw["branchdc_ne"][j]["cost"]))
                revenue=revenue+result_mip["solution"]["nw"][i]["branchdc_ne"][j]["pt"]*(nw["gen"][string(b["tbusdc"])]["cost"][1]-nw["gen"][string(b["fbusdc"])]["cost"][1])*mn_data_mip["scenario_prob"]["1"]
                #capex_cable=capex_cable+nw["branchdc_ne"][j]["cost"]*mn_data_mip["scenario_prob"]["1"]
            end
        end
        for (j,c) in nw["convdc"]

            #capex_conv=capex_conv+c["cost"]*result_mip["solution"]["nw"][i]["convdc"][j]["p_pacmax"]*mn_data_mip["scenario_prob"]["1"]/3
            #capex_conv=capex_conv+c["cost"]
        end
        #wnd_energy=wnd_energy+result_mip["solution"]["nw"][i]["gen"]["2"]["pg"]*mn_data_mip["scenario_prob"]["1"]
        #wnd_revenue_min=wnd_revenue_min+result_mip["solution"]["nw"][i]["gen"]["2"]["pg"]*min(nw["gen"]["1"]["cost"][1])*mn_data_mip["scenario_prob"]["1"]
        wnd_energy=wnd_energy+result_mip["solution"]["nw"][i]["gen"]["4"]["pg"]*mn_data_mip["scenario_prob"]["1"]
        wnd_revenue_min=wnd_revenue_min+result_mip["solution"]["nw"][i]["gen"]["4"]["pg"]*min(nw["gen"]["1"]["cost"][1],nw["gen"]["2"]["cost"][1],nw["gen"]["3"]["cost"][1])*mn_data_mip["scenario_prob"]["1"]
        wnd_revenue_max=wnd_revenue_max+result_mip["solution"]["nw"][i]["gen"]["4"]["pg"]*max(nw["gen"]["1"]["cost"][1],nw["gen"]["2"]["cost"][1],nw["gen"]["3"]["cost"][1])*mn_data_mip["scenario_prob"]["1"]
        wnd_revenue_hm=wnd_revenue_hm+result_mip["solution"]["nw"][i]["gen"]["4"]["pg"]*nw["gen"]["2"]["cost"][1]*mn_data_mip["scenario_prob"]["1"]
    end
    Y_=5/30
    println(string(revenue*Y_)*" "*
    string(capex_cable)*" "*
    string(capex_conv)*" "*
    string(revenue+capex_cable+capex_conv)*" "*
    string(-1*result_mip["objective"])*" "*
    string(wnd_revenue_min*Y_)*" "*
    string(wnd_revenue_hm*Y_)*" "*
    string(wnd_revenue_max*Y_)*" "*
    string(wnd_energy*(100/(length(scenario_years)*624))*8760*scenario_planning_horizon*Y_))
