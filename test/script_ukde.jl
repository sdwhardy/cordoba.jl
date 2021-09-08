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

ic_mva=500;#MVA - interconnector
owpp_mva=1000#mva
ic_length=600;#km
owpp_km=300#km
candidates=[[0.5,0.4],[0.5,0.4]]#[[ic],[owpp]]=#
d1="BEdata";d2="UKdata"
#800MWh with ic_mva=500;owpp_mva=1000#mva,ic_length=600;#km,owpp_km=300#km,candidates=[[0.5,0.4],[0.5,0.4]],x=3 in storage_costs
##############################
function cordoba_go!(d1,d2,ic_mva,owpp_mva,ic_length,owpp_km,candidates::Vector{Array{Float64,1}}=[[1.0],[1.0]])
    casename = "ic"
    file = "./test/data/input/$casename.m"
        data = _PM.parse_file(file)#load data in PM format
        data=_CBD.additional_candidates(data,candidates)
        _PMACDC.process_additional_data!(data)#add extra DC model data
        ##################### Network cost/elec PARAMETERS ######################
        hoa=_CBD.hoa_datastruct_candidateICs(ic_mva,owpp_mva,ic_length,owpp_km,candidates)
        data=_CBD.pm_dict_candidateICs(data,hoa,candidates)
        converter_parameters_rxb(data)#adjust converter parameters
        ##########################################################################

        #################### Multi-period INPUT PARAMETERS #######################
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
        data, z0,z1 = _CBD.get_profile_data_sets(d1,d2,data, n, scenario)
        n=number_of_hours = length(z0.time_stamp)
        #set problem dimension
        dim = scenario["hours"] * length(data["scenario"])

         # create a dictionary to pass time series data to data dictionary
        extradata = _CBD.create_profile_sets(dim, data, z0,z1,ic_mva,owpp_mva)

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
        resultACDC = _PMACDC.run_mp_acdctnepopf(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

    return resultACDC, mn_data
end
