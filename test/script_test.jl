################## loads external packages ##############################
using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV
import cordoba; const _CBD = cordoba#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import FlexPlan; const _FP = FlexPlan
import InfrastructureModels; const _IM = InfrastructureModels
include("../test/data/conv_spec.jl")

################### ENTSO-E scenario ####################################
_sc=["0","3","6"]
_yr="2020"
################### Input data for script_test.jl ########################
owpp_mva=4000#mva of wf (in de)
#interconnectors format: (mva,(km,type)) type: -1:both ends onshore, 0:both ends offshore, 1:onshore/offshore
#must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
ics=[(2000,(550,-1)),(2000,(760,-1)),(2000,(250,-1)),(2000,(470,1)),(2000,(145,1)),(2000,(246,1))];
#location of nodes
zs=["UK","DE","DK","DE"]#must be in same order as .m file gens
candidate_ics=[1,4/5,3/5,1/2]#Candidate Cables
infinite_grid=sum(first.(ics))+owpp_mva#ensures enough generation and consumption in all markets

################## reads .m input file name ######################
casename = "test"
file = "./test/data/input/$casename.m"
data = _PM.parse_file(file)#load data in PM format
data=_CBD.additional_candidatesICS(data,candidate_ics,ics)#adds additional candidates

####################### These functions are not needed if data is directly entered into .m file but this is a recipe for disaster
#essentially it is just duplicating connections in branchdc_ne and adding real cable values for power flow
#Cost of a candidate is the converter/plarforms/transformers and cables all lumped together
#################### Calculates cables for DC lines
for (i,bdc) in data["branchdc_ne"]
data["branchdc_ne"][i]=_CBD.candidateIC_cost(bdc);end
data["branchdc_ne"]=_CBD.unique_candidateIC(data["branchdc_ne"])#keep only unique candidates

iran=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]
iran=["1","2","4","5","7","8","10","11","12","13","14","16","18","19","20","21","22","23","24"]


for br in iran
delete!(data["branchdc_ne"],br);end

for (key,bdc) in data["branchdc_ne"]; bdc["cost"]=0; end
################################################################################

#################### PowerModelsACDC settings ############################
_PMACDC.process_additional_data!(data)#add extra DC model data
converter_parameters_rxb(data)#adjust converter parameters

#################### Multi-period input parameters #######################
n=number_of_hours = 8760 # Number of time points in DE
scenario = Dict{String, Any}("hours" => number_of_hours, "sc_years" => Dict{String, Any}())
scenario["sc_years"]["1"] = Dict{String, Any}()
scenario["sc_years"]["1"]["year"] = 2020#year of data
scenario["sc_years"]["1"]["start"] = Dates.datetime2epochms(DateTime(scenario["sc_years"]["1"]["year"]))   #seconds from year 0 (2020 = 63745056000000)
scenario["sc_years"]["1"]["probability"] = 1
scenario["planning_horizon"] = 25 # in years, to scale generation cost

##################### load time series data ##############################
#time series data files are UKdata.csv, DEdata.csv, DKdata.csv
#scenario
data,zs_data =  _CBD.combine_profile_data_sets_entso_scenario(unique(zs),data, n,_sc,_yr, scenario)

#single set
#data, zs_data = _CBD.get_profile_data_sets_entso_scenario(unique(zs),data, n,_sc[3],_yr, scenario); unique!(zs_data,:time_stamp)#

#sample
zs_data=_CBD.set_of_hours(zs_data,[i for i=1:1:24],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31],[1,2,3,4,5,6,7,8,9,10,11,12])
#########################################################################


#set problem dimension
n=number_of_hours =scenario["hours"]= length(zs_data.time_stamp)
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

#optimation settings
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)

#run optimization
resultACDC = _CBD.cordoba_mp_acdctnepopf(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)


profit_1to2=0
profit_1to3=0
profit_2to3=0
profit_1to4=0
profit_2to4=0
profit_3to4=0
wnd_energy=0
wnd_profit=0
capex=0
mn_data["nw"]["1"]["gen"]
for (i,n) in resultACDC["solution"]["nw"]
    #if (parse(Int,i)>5*8360 && parse(Int,i)<6*8360)
        wnd_energy=wnd_energy+n["gen"]["4"]["pg"]
        wnd_profit=wnd_profit+n["gen"]["4"]["pg"]*max(mn_data["nw"][i]["gen"]["1"]["cost"][1],mn_data["nw"][i]["gen"]["2"]["cost"][1],mn_data["nw"][i]["gen"]["3"]["cost"][1])
        profit_1to2=profit_1to2+n["branchdc_ne"]["3"]["pt"]*(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["1"]["cost"][1])
        profit_1to3=profit_1to3+n["branchdc_ne"]["6"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["1"]["cost"][1])
        profit_2to3=profit_2to3+n["branchdc_ne"]["9"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["2"]["cost"][1])
        profit_1to4=profit_1to4+n["branchdc_ne"]["15"]["pt"]*(mn_data["nw"][i]["gen"]["1"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        profit_2to4=profit_2to4+n["branchdc_ne"]["17"]["pt"]*(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        #profit_3to4=profit_3to4+n["branchdc_ne"]["23"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        capex=capex+(mn_data["nw"][i]["branchdc_ne"]["3"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["6"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["9"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["15"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["17"]["cost"])#+mn_data["nw"][i]["branchdc_ne"]["23"]["cost"])
    #end
end
println(string(profit_1to2)*" "*
string(profit_1to3)*" "*
string(profit_2to3)*" "*
string(profit_1to4)*" "*
string(profit_2to4)*" "*
string(profit_3to4)*" "*
string(capex)*" "*
string(resultACDC["objective"])*" "*
string(wnd_energy)*" "*
string(wnd_profit)*" "*
string(8360)*" "*
string(scenario["planning_horizon"]))




uk2dkwind=deepcopy(resultACDC)
uk2dkwind_dta=deepcopy(mn_data)


for (key,nw) in resultACDC["solution"]["nw"]
    for (key2,b) in nw["ne_storage"]
    #println(key*" -  "*key2*" dr:  "*string(b["sd_ne"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["discharge_rating"]))
    #println(key*" -  "*key2*" cr: "*string(b["sc_ne"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["charge_rating"]))
    #println(key*" -  "*key2*" cr: "*string(b["sc_ne"]*mn_data["nw"][key]["ne_storage"][key2]["cost_abs"])*" "*" dr:  "*string(b["sd_ne"]*mn_data["nw"][key]["ne_storage"][key2]["cost_inj"]))
    println(key*" -  "*key2*" cst: "*string(mn_data["nw"][key]["ne_storage"][key2]["inst_cost"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["eq_cost"]))
    #println(key*" -  "*key2*" cst: "*string(mn_data["nw"][key]["ne_storage"][key2]["cost_abs"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["cost_inj"]))
end;end


using Plots
plotly()
width=750
height=500
p=plot(size = (width, height),xaxis = ("Date", font(20, "Courier")),yaxis = ("P.U.", font(20, "Courier")))
plot!(p,zs_data2020["time_stamp"],zs_data2020["Wnd_MWhDE"],color = :red, label="",size = (width, height))
plot!(p,zs_data2019["time_stamp"],zs_data2019["Wnd_MWhDE"],color = :blue, label="",size = (width, height))

gui()
zs_data2020["Wnd_MWhDE"]=deepcopy(zs_data)
zs_data2019=deepcopy(zs_data)
