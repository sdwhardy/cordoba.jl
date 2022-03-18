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
    scenario_names=["EU17","EU18","EU19","EU20"]
    #scenario_names=["EU19"]

    #scenario_years=["2020","2030","2040"]
    #scenario_years=["2020","2030"]
    scenario_years=["2020"]

    scenario_planning_horizon=30

    ################### Input data for script_test.jl ########################
    owpp_mva=[4000]#mva of wf (in de)
    #owpp_mva=[0]#mva of wf (in de)

    #interconnectors format: (mva,km)
    ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];
    #ics=[(2000,550),(2000,760),(2000,250),(0,470),(0,145),(0,246)];
    #ics=[(4000,145)];
    conv_lim=4000

    #location of nodes
    markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
    #markets_wfs=[["DE"],["DE"]]
    infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

    ##################### load time series data ##############################
    k=2
    #=scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series

    ##################### Cluster time series data ###########################

    for (sc,yrs_ts) in scenario_data
        for (yr,ts) in yrs_ts
            println(sc*" "*yr)

            filter!(row -> ismissing(row.time_stamp)==false, scenario_data[sc][yr])
            filter!(row -> ismissing(row.Wnd_MWhDE)==false, scenario_data[sc][yr])
            filter!(row -> ismissing(row.EUR_daDE)==false, scenario_data[sc][yr])
            filter!(row -> ismissing(row.EUR_daUK)==false, scenario_data[sc][yr])
            filter!(row -> ismissing(row.EUR_daDK)==false, scenario_data[sc][yr])
            daily_ts=_CBD.daily_tss(DateTime.(ts[!,"time_stamp"]))
            daily_dew=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhDE"]))
            daily_de=_CBD.daily_tss(Float64.(ts[!,"EUR_daDE"]))
            daily_uk=_CBD.daily_tss(Float64.(ts[!,"EUR_daUK"]))
            daily_dk=_CBD.daily_tss(Float64.(ts[!,"EUR_daDK"]))

            kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2).+(daily_dk.^2).+(daily_dew.^2))',axis=1), k)
            ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
            filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
        end
    end
    save("./test/data/input/EUSTDG17t020_TS_k"*string(k)*".jld2",scenario_data)=#
    scenario_data=load("./test/data/input/EUSTDG17t020_TS_k"*string(k)*".jld2")
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],scenario_years));delete!(scenario_data[k],y);end; end;end;end
    #d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);end;end
    #scenario_data["EU19"]["2020"]=scenario_data["EU19"]["2020"][1:2,:]

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
    #println(scenario_data)
    #save("./test/data/input/EUSTDG_TS_k"*string(k)*".jld2",scenario_data)

    #scenario_data["EU19"]["2020"]=scenario_data["EU19"]["2020"][1:2,:]



    ################################################# MIP TNEP #######################################
    ################## reads .m input file name ######################
    #casename = "test_convex_wf"
    #casename = "test_convset_wf"
    #casename = "test_convex_conv"
    casename = "test_convex_conv_set_nippon_convexcable"
    file = "./test/data/input/$casename.m"
    data_mip = _PM.parse_file(file)#load data in PM format

    #################### Calculates cables for DC lines
    z_base_dc=(data_mip["busdc"]["1"]["basekVdc"])^2/data_mip["baseMVA"]
    genz=[];wfz=[]
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid/data_mip["baseMVA"]));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid/data_mip["baseMVA"]));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]/data_mip["baseMVA"]));end
    #candidate_ics=[1]#Candidate Cable sizes
    #candidate_ics=[4/5]#Candidate Cable sizes
    candidate_ics=[1,4/5,3/5,1/2,1/4,1/5]#Candidate Cable sizes
    #candidate_ics=[1/2]#Candidate Cable sizes
    #candidate_ics=[1]#Candidate Cable sizes
    data_mip=_CBD.additional_candidatesICS_DC(data_mip,candidate_ics,ics)#adds additional candidates
    for (i,bdc) in data_mip["branchdc_ne"]
    data_mip["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance_DC(bdc,z_base_dc);end
    #data_mip["branchdc_ne"]=_CBD.unique_candidateIC_DC(data_mip["branchdc_ne"])#keep only unique candidates

    ######### delete rename cables - set up the problem
    #data_mip["branchdc_ne"]["13"]["rateA"]=0;data_mip["branchdc_ne"]["14"]["rateA"]=0;;data_mip["branchdc_ne"]["15"]["rateA"]=0;data_mip["branchdc_ne"]["16"]["rateA"]=0
    #delete!(data_mip["branchdc_ne"],"4");delete!(data_mip["branchdc_ne"],"5");delete!(data_mip["branchdc_ne"],"6")
    #=delete!(data_mip["branchdc_ne"],"7");delete!(data_mip["branchdc_ne"],"8");delete!(data_mip["branchdc_ne"],"9")
    delete!(data_mip["branchdc_ne"],"10");delete!(data_mip["branchdc_ne"],"11");delete!(data_mip["branchdc_ne"],"12")
    delete!(data_mip["branchdc_ne"],"14");delete!(data_mip["branchdc_ne"],"15");delete!(data_mip["branchdc_ne"],"16")
    delete!(data_mip["branchdc_ne"],"19");delete!(data_mip["branchdc_ne"],"18");delete!(data_mip["branchdc_ne"],"20")
    delete!(data_mip["branchdc_ne"],"23");delete!(data_mip["branchdc_ne"],"22");delete!(data_mip["branchdc_ne"],"24")=#
    _PMACDC.process_additional_data!(data_mip)#add extra DC model data
    _CBD.converter_parameters_rxb(data_mip)

    #################### Multi-period input parameters #######################
    all_scenario_data,data_mip,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_mip)
    scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
    extradata,data_mip = _CBD.create_profile_sets_mesh(dim, data_mip, all_scenario_data, markets_wfs, infinite_grid, data_mip["baseMVA"])

    #################### Scale cost data
    #_CBD.scale_cost_data_2hourly!(data_mip, scenario)
    #_CBD.scale_cost_data_2hourly!(extradata, scenario)
    _CBD.scale_cost_data_2hourly!(data_mip, scenario)
    _CBD.scale_cost_data_2yearly!(extradata, scenario)

    # Create data dictionary where time series data is included at the right place
    mn_data_mip = _PMACDC.multinetwork_data(data_mip, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()
    mn_data_mip= _CBD.npvs_costs_datas_4mip(mn_data_mip, scenario, scenario_years)
    #[println(k*" "*b*" "*string(br["cost"])) for (k,nw) in mn_data_mip["nw"] for (b,br) in nw["convdc"]];println()

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
#-14571.410431569195 16,20,4
    admm_4_AjAwAgAuAo_main(mn_data_mip, gurobi, s)

    #result_mip = _CBD.cordoba_mp_acdctnepopf_stoch(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    #result_mip = _CBD.acdc_tnep_convex_conv_npv(mn_data_mip, _PM.LPACCPowerModel, gurobi, multinetwork=true; setting = s)
    #result_mip = _CBD.acdc_tnep_convex_conv_strg_npv(mn_data_mip, _PM.SOCWRPowerModel, gurobi, multinetwork=true; setting = s)

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
s=deepcopy(s2)
r2=deepcopy(results_set)

#verify
a="Aox"
s["agent"]=a
result_mip = _CBD.acdc_tnep_convex_conv_strg_admm_convexcable(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
s["fixed_variables"] = _CBD.admm_4_AjAwAgAuAo_update(result_mip["solution"]["nw"], s["fixed_variables"], a, wfz)
a="Aobx"
s["agent"]=a
result_mip = _CBD.acdc_tnep_convex_conv_strg_admm_convexcable(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#try to print objective function
s["fixed_variables"] = _CBD.admm_4_AjAwAgAuAo_update(result_mip["solution"]["nw"], s["fixed_variables"], a, wfz)
