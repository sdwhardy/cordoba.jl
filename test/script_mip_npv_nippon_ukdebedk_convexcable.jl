    ################## loads external packages ##############################
    using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, Clp
    using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    import FlexPlan; const _FP = FlexPlan
    #import InfrastructureModels; const _IM = InfrastructureModels
    include("../aux/post_process/functions.jl")

    ################### ENTSO-E scenario description ####################################
    #scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
    #scenario_names=["EU17","EU18","EU19","EU20"]
    scenario_names=["EU19"]

    #scenario_years=["2020","2030","2040"]
    #scenario_years=["2020","2030"]
    scenario_years=["2030"]

    scenario_planning_horizon=30

    ################### Input data for script_test.jl ########################
    #owpp_mva=[4000]
    owpp_mva=[4000,4000]#mva of wf (in de)
    #owpp_mva=[0]#mva of wf (in de)

    #interconnectors format: (mva,km)
    #ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];
    ics=[(4000,550),(4000,760),(4000,250),(4000,205),(4000,596),(4000,647),
    (4000,470),(4000,145),(4000,246),(4000,160),(4000,45),(4000,602),(4000,462)];
    #ics=[(2000,550),(2000,760),(2000,250),(0,470),(0,145),(0,246)];
    #ics=[(4000,145)];
    conv_lim=4000

    #location of nodes
    #markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
    markets_wfs=[["UK","DE","DK","UK","BE"],["DE","BE"]]#must be in same order as .m file gens
    #markets_wfs=[["DE"],["DE"]]
    infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

    ##################### load time series data ##############################
    k=2
    scenario_data=load("./test/data/input/BE_DE_UK_DK_islands/EUSTDG_17to20_TS_k"*string(k)*"_wBE.jld2")
    #take only desired scenarios and years
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],scenario_years));delete!(scenario_data[k],y);end; end;end;end

    ##################### Find minimum length scenario
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    ##################### Make all scenarios the same length
    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end



    ################################################# MIP TNEP #######################################
    ################## reads .m input file name ######################
    casename = "test_convex_conv_set_nippon_ukdebedk"
    file = "./test/data/input/$casename.m"
    data_mip = _PM.parse_file(file)#load data in PM format

    #################### Calculates cables for DC lines
    z_base_dc=(data_mip["busdc"]["1"]["basekVdc"])^2/data_mip["baseMVA"]
    genz=[];wfz=[]#store conventional gens and wind seperately
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid/data_mip["baseMVA"]));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid/data_mip["baseMVA"]));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]/data_mip["baseMVA"]));end

    #Candidate Cable sizes
    candidate_ics=[1,4/5,3/5,1/2]#Candidate Cable sizes

    data_mip=_CBD.additional_candidatesICS(data_mip,candidate_ics,ics)#adds additional candidates
    for (i,bdc) in data_mip["branchdc_ne"]
    data_mip["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance(bdc,z_base_dc);end
    data_mip["branchdc_ne"]=_CBD.unique_candidateIC(data_mip["branchdc_ne"])#keep only unique candidates

    ######### delete/rename cables - set up adjusted problem
    #data_mip["branchdc_ne"]["13"]["rateA"]=0;data_mip["branchdc_ne"]["14"]["rateA"]=0;;data_mip["branchdc_ne"]["15"]["rateA"]=0;data_mip["branchdc_ne"]["16"]["rateA"]=0
    #delete!(data_mip["branchdc_ne"],"4");delete!(data_mip["branchdc_ne"],"5");delete!(data_mip["branchdc_ne"],"6")
    #=delete!(data_mip["branchdc_ne"],"7");delete!(data_mip["branchdc_ne"],"8");delete!(data_mip["branchdc_ne"],"9")
    delete!(data_mip["branchdc_ne"],"10");delete!(data_mip["branchdc_ne"],"11");delete!(data_mip["branchdc_ne"],"12")
    delete!(data_mip["branchdc_ne"],"14");delete!(data_mip["branchdc_ne"],"15");delete!(data_mip["branchdc_ne"],"16")
    delete!(data_mip["branchdc_ne"],"19");delete!(data_mip["branchdc_ne"],"18");delete!(data_mip["branchdc_ne"],"20")
    delete!(data_mip["branchdc_ne"],"23");delete!(data_mip["branchdc_ne"],"22");delete!(data_mip["branchdc_ne"],"24")=#

    #print_topology_data(data_mip)
    _PMACDC.process_additional_data!(data_mip)#add extra DC model data
    _CBD.converter_parameters_rxb(data_mip)

    #################### Multi-period input parameters #######################
    all_scenario_data,data_mip,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_mip)
    scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
    extradata,data_mip = _CBD.create_profile_sets_mesh(dim, data_mip, all_scenario_data, markets_wfs, infinite_grid, [data_mip["baseMVA"],data_mip["baseMVA"]])
    print_topology_data(data_mip,markets_wfs)

    #################### Scale cost data (no NPV) accounted for
    _CBD.scale_cost_data_2hourly!(data_mip, scenario)
    _CBD.scale_cost_data_2yearlyhourly!(extradata, scenario)

    # Create data dictionary where time series data is included at the right place
    mn_data_mip = _PMACDC.multinetwork_data(data_mip, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))

    #scale data considering NPV (4% assumed)
    mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)
    mn_data_mip= _CBD.npvs_costs_datas_4mip(mn_data_mip, scenario, scenario_years)

    #run optimization
    #optimation settings
    #, "home_market" => [2,4]
    s = Dict("output" => Dict("branch_flows" => false, "duals"=>false),"fixed_variables" => Dict{String,Any}(),"agent" => "","relax_problem" => false,
    "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls,
    "corridor_limit" => true, "genz"=>genz, "wfz"=>wfz, "ic_lim"=>conv_lim/100, "strg_lim_offshore"=>0.2, "strg_lim_onshore"=>10, "compare_mode" => false, "dual_update"=>false)
    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
    #
    result_mip = _CBD.acdc_tnep_convex_conv_strg_npv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)# -22240.628145837014 24-12-20-13, -21908.740816338737 4-5-6
    print_solution_data(result_mip, data_mip)
    admm_4_AjAwAgAuAo_main(mn_data_mip, gurobi, s)

    #result_mip = _CBD.cordoba_mp_acdctnepopf_stoch(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    #result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.LPACCPowerModel, gurobi, multinetwork=true; setting = s)
    #result_mip = _CBD.acdc_tnep_convex_conv_strg_npv(mn_data_mip, _PM.SOCWRPowerModel, gurobi, multinetwork=true; setting = s)
#-68673.52170154548 9,10,13,4, 11, 7
s["fixed_variables"]["1"]
a=agents[1]
a=agents[2]
a=agents[3]
a=agents[4]
###########################  development ########################
function admm_4_AjAwAgAuAo_main(mn_data_mip, gurobi, s)
    results_set=[]
    eps=1#set max value of residual for convergence
    residual=Inf
    agents=["Au","Ag","Aw","Aj","Ao"]#,"Ag","Au","Aw"]#,"Aox"]#,"Aw","Aj"]#Initilize agants
    s["fixed_variables"] = _CBD.admm_4_AjAwAgAuAo_intialize(mn_data_mip["nw"], s["fixed_variables"], s["wfz"], s["genz"])
    while (residual>eps)
    #for i=1:10
        results_set=[]
        for a in agents
            s["agent"]=a
            println(a)
            result_mip = _CBD.acdc_tnep_convex_conv_strg_admm_convexcable(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
            s["fixed_variables"] = _CBD.admm_4_AjAwAgAuAo_update(result_mip["solution"]["nw"], s["fixed_variables"], a, wfz)
            push!(results_set,(a,result_mip))
        end
        [println(first(r)*" "*string(last(r)["objective"])) for r in results_set]
        s["fixed_variables"], residual = _CBD.dual_variable_update(s["fixed_variables"])
        println("Residual: "*string(residual))
    #end
    end
    return result_mip
end

#verify
a="Aox"
s["agent"]=a
result_mip = _CBD.acdc_tnep_convex_conv_strg_admm_convexcable(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
s["fixed_variables"] = _CBD.admm_4_AjAwAgAuAo_update(result_mip["solution"]["nw"], s["fixed_variables"], a, wfz)
a="Aobx"
s["agent"]=a
result_mip = _CBD.acdc_tnep_convex_conv_strg_admm_convexcable(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
for (b,br) in result_mip["solution"]["nw"]["1"]["branchdc_ne"];if (br["isbuilt"]==1);println(b*" "*string(br["isbuilt"]));end;end

#=
#-63942.25484022916
10 1.0
44 1.0
20 1.0
6 1.0
4 1.0
16 1.0
36 1.0
48 1.0

-"objective" => -7.04e+4…
32 1.0
42 1.0
28 1.0
37 1.0
12 1.0
16 1.0
36 1.0
3 1.0=#
