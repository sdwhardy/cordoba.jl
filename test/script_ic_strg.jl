using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV
import cordoba; const _CDB = cordoba
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import InfrastructureModels; const _IM = InfrastructureModels
include("../src/economics/main.jl")
include("../src/cordoba/cordoba.jl")
include("../data/conv_spec.jl")


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE EENS is set to zero for flexlan !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#############################
#=data["convdc_ne"]["1"]["cost"]=0
data["convdc_ne"]["2"]["cost"]=0
data["convdc_ne"]["3"]["cost"]=0
data["branchdc_ne"]["1"]["cost"]=0
data["branchdc_ne"]["2"]["cost"]=0
data["branchdc_ne"]["3"]["cost"]=1e10
data["ne_branch"]["1"]["construction_cost"]=0
data["ne_branch"]["2"]["construction_cost"]=0
mn_data["nw"]["1"]["gen"]["1"]["cost"]
mn_data["nw"]["1"]["gen"]["2"]["cost"]
for (i,nw) in mn_data["nw"]
    if (nw["gen"]["4"]["cost"][1] != nw["gen"]["1"]["cost"][1])
    println(i*": "*string(nw["gen"]["4"]["cost"])*" - "*string(nw["gen"]["1"]["cost"]));end
end
mn_data["nw"]["31"]["gen"]["4"]["cost"]=#
ic_mva=500;#MVA - interconnector
owpp_mva=1000#mva
ic_length=600;#km
owpp_km=300#km
candidates=[[0.5,0.4],[0.5,0.4]]#[[ic],[owpp]]=#
d1="BEdata";d2="UKdata"
#800MWh with ic_mva=500;owpp_mva=1000#mva,ic_length=600;#km,owpp_km=300#km,candidates=[[0.5,0.4],[0.5,0.4]],x=3 in storage_costs
##############################
function cordoba_go!(d1,d2,ic_mva,owpp_mva,ic_length,owpp_km,candidates::Vector{Array{Float64,1}}=[[1.0],[1.0]])
    casename = "cordoba_candidateic_wstrg"
    file = "./test/data/input/$casename.m"
    data = _PM.parse_file(file)#load data in PM format
    data=additional_candidates(data,candidates)
    _PMACDC.process_additional_data!(data)#add extra DC model data
    _FP.add_storage_data!(data) # Add addtional storage data model
    ##################### Network cost/elec PARAMETERS ######################
    hoa=hoa_datastruct_candidateICs(ic_mva,owpp_mva,ic_length,owpp_km,candidates)
    data=pm_dict_candidateICs(data,hoa,candidates)
    converter_parameters_rxb(data)#adjust converter parameters
    storage_costs(data,2030)#adjust storage parameters - default: 2021 also 2030 included
    ##########################################################################

    #################### Multi-period INPUT PARAMETERS #######################
    #dn=7;wn=[1,2,3,4,5];mn=12
    #number_of_hours = dn*length(wn)*mn*24 # Number of time points
    n=number_of_hours = 1000 # Number of time points
    scenario = Dict{String, Any}("hours" => number_of_hours, "sc_years" => Dict{String, Any}())
    scenario["sc_years"]["1"] = Dict{String, Any}()
    scenario["sc_years"]["1"]["year"] = 2020#year of data
    scenario["sc_years"]["1"]["start"] = Dates.datetime2epochms(DateTime(scenario["sc_years"]["1"]["year"]))   #seconds from year 0 (2020 = 63745056000000)
    scenario["sc_years"]["1"]["probability"] = 1
    scenario["planning_horizon"] = 30 # in years, to scale generation cost
    ###########################################################################
     # get wind, generation and load time series
    #data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses = get_profile_data_UKBE(data, scenario)
    data, z0,z1 = get_profile_data_sets(d1,d2,data, n, scenario)
    n=number_of_hours = length(z0.time_stamp)
    #set problem dimension
    dim = scenario["hours"] * length(data["scenario"])

     # create a dictionary to pass time series data to data dictionary
    #extradata = create_profile_UKBE_0835_update(dim, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses, ic_mva, owpp_mva)
    extradata = create_profile_sets(dim, data, z0,z1)
    #extradata =Dict{String,Any}();push!(extradata,"dim"=>2016)
    # Scale cost data
    _FP.scale_cost_data_cordoba!(data, scenario)
    _FP.scale_cost_data_cordoba!(extradata, scenario)
#35.4375 18.71 0.57 -8.11
    # Scale battery parameters
    #println((scenario["hours"] / (8760*scenario["planning_horizon"]))*data["ne_storage"]["1"]["max_energy_absorption"])
    _FP.scale_bat_data_cordoba!(data, scenario)

    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type","scenario","name", "source_version", "per_unit"]))

    #print values to REPL
    #printEQuipValues(mn_data["nw"]["1"])#7542 is full power at OWPP and IC
    println(number_of_hours)
    [println(string(i)*", er: "*string(s["energy_rating"])*", mea: "*string(s["max_energy_absorption"])*", dr: "*string(s["discharge_rating"])*", cr: "*string(s["charge_rating"])*", ec: "*string(s["eq_cost"])*", ic: "*string(s["inst_cost"])) for (i,s) in mn_data["nw"]["1"]["ne_storage"]]

    #print values to REPL
    #printEQuipValues(mn_data["nw"]["1"])#7542 is full power at OWPP and IC

    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

    #settings
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)

    #run optimization
    #resultDC = run_tnepopf(data, DCPPowerModel, gurobi, setting = s)
    resultACDC = _FP.cordoba_strg_tnep(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    [println(i) for (i,res) in resultACDC["solution"]["nw"]["1"]["ne_storage"] if (res["isbuilt"]==1)]
    println();println()
    [println(i) for (i,res) in resultACDC["solution"]["nw"]["1"]["ne_storage"] if (res["isbuilt"]==1)]

    return resultACDC, mn_data
end
println(z0)
println(z1)
for (k,d) in sort(extradata["ne_storage"],rev=true)
    println(k*": ")
    println("dcr: "*string(d["discharge_rating"]*data["baseMVA"]))
    println("dcc: "*string(d["cost_inj"]*1000000/data["baseMVA"]))
    println("cr: "*string(d["charge_rating"]*data["baseMVA"]))
    println("cc: "*string(d["cost_abs"]*1000000/data["baseMVA"]))
end
#=resultACDC, mn_data= owpp_ic_beuk(ic_mva,owpp_mva,ic_length,owpp_km)
println(resultACDC["objective"])
for (i,g) in resultACDC["solution"]["nw"]["1"]["gen"]
    println(string(i)*" - "*string(g["pg"]))
end=#
#data arrays
#connections=Vector{Int64}();costs=Vector{Float64}();times=Vector{Float64}();ic_mvas=Vector{Float64}();
#ic_kms=Vector{Float64}();owpp_mvas=Vector{Float64}();owpp_kms=Vector{Float64}()
resultACDC["solution"]["nw"]["1856"]
mn_data["nw"]["1856"]=#
for (k,nw) in mn_data["nw"]
    println(k*" - "*string(nw["gen"]["5"]["pmax"]))
end
file="../FlexPlan/test/data/cordoba/results/owpp_ic_nostg_2GW600km_candidateIC.csv"
wcfile=open(file,"w")
println(wcfile,"connection,ic_mva,ic_km,owpp_mva,owpp_km,strg,ic,ac,cnv")
close(wcfile)
for ic_mva in [500]
    #for owpp_mva in [300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500]
    for owpp_mva in [3000]
    #for owpp_mva in [750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500]
    #for owpp_mva in [100]
        for ic_length in [600]
            owpp_kms=[i for i=20:10:ic_length-20]
            #owpp_kms=[i for i=20:10:ic_length-20]
            #owpp_kms=[5]
            for owpp_km in owpp_kms
            println("ic_mva: "*string(ic_mva)," owpp_mva: "*string(owpp_mva)," ic_length: "*string(ic_length)," owpp_km: "*string(owpp_km))
            resultACDC, mn_data=owpp_ic_beuk(ic_mva,owpp_mva,ic_length,owpp_km,[[1.0],[1.0]])#[[ic],[owpp]]
            wcfile=open(file,"a")
            store_data(wcfile,resultACDC, ic_mva, owpp_mva, ic_length, owpp_km)
            close(wcfile)
            println(resultACDC["primal_status"])
            println()
end;end;end;end


#connection (-1:BE, 0:UK, 1:HOA), cost, ic_mva, owpp_mva
#df = DataFrame(connection = connections, cost = costs, time = times, ic_mva=ic_mvas, ic_km=ic_kms, owpp_mva=owpp_mvas, owpp_km=owpp_kms)
#CSV.write("../FlexPlan/test/data/owpp_ic_2be_nostg.csv", df)




################################################################################
#post process
cv, DCbr, ACbr=display_results(resultACDC)
resultACDC["solution"]["nw"]["8350"]
mn_data["nw"]["8350"]["gen"]
resultACDC["solution"]["nw"]["1"]
rc=deepcopy(resultACDC["solution"]["nw"]["1"])
#findmax(windgenprofile_beuk[1,:])
g="2"
mn_data["nw"]["8350"]["gen"][g]["pmax"]
mn_data["nw"]["8350"]["gen"][g]["cost"]












##########################################

for (i,nw) in mn_data["nw"]
    for (g,gen) in nw["gen"]
        println(i*" - "*g*" pmax "*string(gen["pmax"])*" pmin "*string(gen["pmin"])*" cost "*string(gen["cost"][1]))
    end
end

######################### equipment models ##########################
mva=1500.0;km=100;m=17.1
ac_cable=AC_cbl(mva, km)#mva, km - 220kV
dc_cable=DC_cbl(mva, km)#mva, km - 300kV
ac_platform=AC_plat(mva, m)#mva, depth in m (optional)
dc_platform=DC_plat(mva, m)#mva, depth in m (optional)
off_transformer=off_xfm(mva)#mva
on_transformer=on_xfm(mva)#mva
off_converter=off_conv(mva)#mva
on_converter=on_conv(mva)#mva
#println("DC: "*string(dc_cable.costs.cpx_p+dc_cable.costs.cpx_i+on_converter.costs.cpx+off_converter.costs.cpx+on_transformer.costs.cpx_p+on_transformer.costs.cpx_i+off_transformer.costs.cpx_p+off_transformer.costs.cpx_i+dc_platform.costs.cpx))
#println("AC: "*string(ac_cable.costs.cpx_p+ac_cable.costs.cpx_i+on_transformer.costs.cpx_p+on_transformer.costs.cpx_i+off_transformer.costs.cpx_p+off_transformer.costs.cpx_i+ac_platform.costs.cpx))
####################################### Addition

println("Generators: ")
for (k,nb) in data["gen"]
    println(k*" pmax - "*string(nb["pmax"])*" gen_bus- "*string(nb["gen_bus"]))
end
println("AC lines: ")
for (k,nb) in data["ne_branch"]
    println(k*" rate_a- "*string(nb["rate_a"])*" f_bus- "*string(nb["f_bus"])*" to_bus- "*string(nb["t_bus"])*" cost- "*string(nb["construction_cost"]))
end
println("DC lines: ")
for (k,nb) in data["branchdc_ne"]
    println(k*" rateA - "*string(nb["rateA"])*" f_bus- "*string(nb["fbusdc"])*" to_bus- "*string(nb["tbusdc"])*" cost- "*string(nb["cost"]))
end
println("DC candidate converters: ")
for (k,nb) in data["convdc_ne"]
    println(k*" Pacmax - "*string(nb["Pacmax"])*" Pacmin- "*string(nb["Pacmin"])*" busdc_i- "*string(nb["busdc_i"])*" busac_i- "*string(nb["busac_i"])*" cost- "*string(nb["cost"]))
end
println("DC converters: ")
for (k,nb) in data["convdc"]
    println(k*" Pacmax - "*string(nb["Pacmax"])*" Pacmin- "*string(nb["Pacmin"])*" busdc_i- "*string(nb["busdc_i"])*" busac_i- "*string(nb["busac_i"]))
end
