using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV
import cordoba; const _CBD = cordoba
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import FlexPlan; const _FP = FlexPlan
import InfrastructureModels; const _IM = InfrastructureModels
include("../test/data/conv_spec.jl")


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE EENS is set to zero for CORDOBA/FlexPlan!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


owpp_mva=4000#mva additional 2.5GW capacity nearby

#-1:both ends onshore, 0:both ends offshore, 1:onshore/offshore
#ics=[(2000,(615,-1)),(2000,(677,-1)),(2000,(215,-1)),(2000,(480,1)),(2000,(176,1)),(2000,(224,1))];#must be in same order as zs uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3)
#zs=["UK","DE","DK","DE"]#must be in same order as .m file gens

ics=[(2000,(550,-1)),(2000,(760,-1)),(2000,(250,-1)),(2000,(470,1)),(2000,(145,1)),(2000,(246,1)),(4000,(145,1))];#must be in same order as zs uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3)

zs=["UK","DE","DK","DE","DE"]#must be in same order as .m file gens
infinite_grid=sum(first.(ics))+owpp_mva#ensures enough generation and consumption in all markets
candidate_ics=[1,4/5,3/5,1/2]#[[ic=1000,1200,1600,2000],[owpp]]#

#800MWh with ic_mva=500;owpp_mva=1000#mva,ic_length=600;#km,owpp_km=300#km,candidates=[[0.5,0.4],[0.5,0.4]],x=3 in storage_costs
##############################
function cordoba_go!(d1,d2,ic_mva,owpp_mva,ic_length,owpp_km,candidates::Vector{Array{Float64,1}}=[[1.0],[1.0]])
    casename = "ic_ukdedk"
    file = "./test/data/input/$casename.m"
    data = _PM.parse_file(file)#load data in PM format
    data=_CBD.additional_candidatesICS(data,candidate_ics,ics)

    #delete!(data["branchdc_ne"],"6");delete!(data["branchdc_ne"],"7");delete!(data["branchdc_ne"],"8")
    #delete!(data["branchdc_ne"],"5");delete!(data["branchdc_ne"],"6")

    for (i,bdc) in data["branchdc_ne"]#onshore to onshore connections
    data["branchdc_ne"][i]=_CBD.candidateIC_cost(bdc);end

    data["branchdc_ne"]=_CBD.unique_candidateIC(data["branchdc_ne"])#keep only unique candidates


    _PMACDC.process_additional_data!(data)#add extra DC model data
    ##################### Network cost/elec PARAMETERS ######################

    converter_parameters_rxb(data)#adjust converter parameters
    ##########################################################################

    #################### Multi-period INPUT PARAMETERS #######################
    n=number_of_hours = 7151 # Number of time points in DE
    scenario = Dict{String, Any}("hours" => number_of_hours, "sc_years" => Dict{String, Any}())
    scenario["sc_years"]["1"] = Dict{String, Any}()
    scenario["sc_years"]["1"]["year"] = 2020#year of data
    scenario["sc_years"]["1"]["start"] = Dates.datetime2epochms(DateTime(scenario["sc_years"]["1"]["year"]))   #seconds from year 0 (2020 = 63745056000000)
    scenario["sc_years"]["1"]["probability"] = 1
    scenario["planning_horizon"] = 30 # in years, to scale generation cost
    ###########################################################################
     # get wind, generation and load time series
    #data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses = get_profile_data_UKBE(data, scenario)
    data, zs_data = _CBD.get_profile_data_sets_mesh(unique(zs),data, n, scenario)
    n=number_of_hours = length(zs_data.time_stamp)
    #set problem dimension
    dim = scenario["hours"] * length(data["scenario"])

     # create a dictionary to pass time series data to data dictionary
    extradata,data = _CBD.create_profile_sets_mesh(dim, data, zs_data, zs, infinite_grid, owpp_mva)
    # Scale cost data
    _FP.scale_cost_data_cordoba!(data, scenario)
    _FP.scale_cost_data_cordoba!(extradata, scenario)

    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type","scenario","name", "source_version", "per_unit"]))

    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

    #settings
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)

    #run optimization
    resultACDC = _CBD.cordoba_mp_acdctnepopf(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

    #resultACDC = _PMACDC.run_tnepopf(data, _PM.DCPPowerModel, gurobi, setting = s)

    return resultACDC, mn_data
end
display(mn_data["gen"])

profit_125=0
profit_522=0
profit_223=0
profit_321=0
wnd_energy=0
wnd_energy_profit=0
for (i,n) in uk2dkFull["solution"]["nw"]
    wnd_energy=wnd_energy+n["gen"]["5"]["pg"]
    wnd_energy_profit=wnd_energy_profit+n["gen"]["5"]["pg"]*max(uk2dkFull_dta["nw"][i]["gen"]["1"]["cost"][1],uk2dkFull_dta["nw"][i]["gen"]["2"]["cost"][1],uk2dkFull_dta["nw"][i]["gen"]["3"]["cost"][1])
    #profit_125=profit_125+abs(n["branchdc_ne"]["11"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["1"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["5"]["cost"][1])
    #profit_522=profit_522+abs(n["branchdc_ne"]["14"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["5"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["2"]["cost"][1])
    #profit_223=profit_223+abs(n["branchdc_ne"]["8"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["2"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["3"]["cost"][1])
    #profit_321=profit_321+abs(n["branchdc_ne"]["4"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["3"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["1"]["cost"][1])

    #profit_uk2dk=profit_uk2dk+abs(n["branchdc_ne"]["6"]["pt"])*abs(mn_data["nw"][i]["gen"]["1"]["cost"][1]-mn_data["nw"][i]["gen"]["3"]["cost"][1])
    #profit_de2dk=profit_de2dk+abs(n["branchdc_ne"]["10"]["pt"])*abs(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["3"]["cost"][1])
    end
uk2dkFull=deepcopy(resultACDC)
uk2dkFull_dta=deepcopy(mn_data)
