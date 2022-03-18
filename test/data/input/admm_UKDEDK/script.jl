    ################## loads external packages ##############################
    #using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames, XLSX
    #using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    using Gurobi, DataFrames, JuMP
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    ########################################################################

    ##################### Input parameters #################################
    rt_ex=pwd()*"\\test\\data\\input\\admm_UKDEDK\\"#folder path
    argz = Dict(
    "scenario_planning_horizon"=>30,
    "scenario_names"=>["EU17","ST18","ST19","DG20"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
    "k"=>2,#number of representative days modelled (24 hours per day)
    "scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
    "owpp_mva"=>[4000],#mva of wf in MVA
    "conv_lim"=>4000,#Max Converter size in MVA
    "candidate_ics_ac"=>[6/5,1,4/5,3/5,2/5,1/5],#AC Candidate Cable sizes (fraction of full MVA)
    "candidate_ics_dc"=>[6/5,1,4/5,3/5,2/5,1/5],#DC Candidate Cable sizes (fraction of full MVA)
    "dr"=>0.04,#discount rate
    "yearly_investment"=>1000000)

    ################## optimization/solver setup options ###################

    s = Dict("output" => Dict("branch_flows" => false),
    #"home_market"=>[2,4],#nodes within HM -25786.84091162341
    #"balancing_reserve"=>0.3,#
    "eps"=>0.1,#admm residual
    "AC"=>"1",#0=false, 1=true
    "relax_problem" => true,
    "conv_losses_mp" => false,
    "process_data_internally" => false,
    "corridor_limit" => true,
    "strg_lim_offshore"=>0.2,
    "strg_lim_onshore"=>10,
    "max_invest_per_year"=>_CBD.max_invest_per_year(argz))

    mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)
    #data["branchdc"]
    #"br_x" => 0.005351380868391557
    #"br_r" => 0.000871…
    #"br_x" => 0.0400
    #"br_r" => 0.00400
    ########################################################################

    ################## optimization/solver setup options ###################

    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)
    result_mip = _CBD.cordoba_acdc_wf_strg_admm(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#-16164.432804419186 21-13-17
    _CBD.print_solution_data(result_mip, data, argz)

    _CBD.print_convex_solution_data(result_mip, data, argz)#-25273.236946867946
    mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)




    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)
    result_mip = _CBD.cordoba_acdc_wf_strg_admm(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#-16164.432804419186 21-13-17
    _CBD.print_solution_data(result_mip, data, argz)

    #ADMM Set up
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
    result_mip = _CBD.admm_4_AjAwAgAuAo_main(mn_data, gurobi, s)
    mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)
    result_mip = _CBD.cordoba_acdc_wf_strg_admm(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#-16164.432804419186 21-13-17
    _CBD.print_solution_data(result_mip, data, argz)





    _CBD.print_solution_data(result_mip, data, argz)
    _CBD.print_convex_solution_data(result_mip, data, argz)
    _CBD.print_convex_solution_data(last(results_set[5]), data, argz)


s["fixed_variables"]["1"]
a=agents[1]
a=agents[2]
a=agents[3]
a=agents[4]
a=agents[5]
###########################  development ########################


########################## gomi
result_mip["solution"]["nw"]["142"]["storage"]
cost_1=0
cost_2=0
cost_3=0
cost_4=0
cost_5=0
cost_6=0
cost_7=0
bat_1_in=0
bat_2_in=0
bat_3_in=0
bat_4_in=0
bat_1_out=0
bat_2_out=0
bat_3_out=0
bat_4_out=0

using OrderedCollections
result_mip["solution"]["nw"]=sort!(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]), by=x->parse(Int64,x))
for (i,nw) in (result_mip["solution"]["nw"])
    #if (parse(Int64,i)>0 && parse(Int64,i)<49)
        cost_1=cost_1+nw["gen"]["1"]["pg"]
        cost_2=cost_2+nw["gen"]["2"]["pg"]
        cost_3=cost_3+nw["gen"]["3"]["pg"]
        cost_4=cost_4+nw["gen"]["4"]["pg"]
        cost_5=cost_5+nw["gen"]["5"]["pg"]
        cost_6=cost_6+nw["gen"]["6"]["pg"]
        cost_7=cost_7+nw["gen"]["7"]["pg"]
        bat_1_in=bat_1_in+nw["storage"]["1"]["ps"]
        bat_2_in=bat_2_in+nw["storage"]["2"]["ps"]
        bat_3_in=bat_3_in+nw["storage"]["3"]["ps"]
        bat_4_in=bat_4_in+nw["storage"]["4"]["ps"]
        println(i)
        println(nw["gen"]["1"]["pg"]+nw["gen"]["2"]["pg"]+nw["gen"]["3"]["pg"]+nw["gen"]["4"]["pg"]+nw["gen"]["5"]["pg"]+nw["gen"]["6"]["pg"]+nw["gen"]["7"]["pg"]-nw["storage"]["1"]["ps"]-nw["storage"]["2"]["ps"]-nw["storage"]["3"]["ps"]-nw["storage"]["4"]["ps"])
        #cost_dc=cost_dc+nw["branchdc"]["1"]["cost"]
    #end
end
println(cost_1)#7444.97039870115 (relax)
println(cost_2)#7050.082286898671 (relax)
println(cost_3)#8660.616524895151 (relax)
println(cost_4)#2345.185999790022 (relax)
println(cost_5)#-10684.822408001537 (relax)
println(cost_6)#-9305.79039924618 (relax)
println(cost_7)#-9596.589778667829 (relax)
println(cost_1+cost_2+cost_3+cost_4+cost_5+cost_6+cost_7-bat_1_in-bat_2_in-bat_3_in-bat_4_in)
