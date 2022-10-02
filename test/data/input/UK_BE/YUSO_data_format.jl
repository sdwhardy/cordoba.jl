using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

#load smaple file
scenario_data_file="C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2"
scenario_data=FileIO.load(scenario_data_file)  
#load YUSO data
df=CSV.read("C://Users//shardy//Documents//julia//times_series_input_large_files//YUSO_data//Cordoba_WP2_data_wTotalCap.csv", DataFrame)
#keep only hourly input
filter!(:dtutc=>x->parse(Int64,x[15:16])==0.0, df);
#trim off sub second values
df.dtutc=[t[1:19] for t in df.dtutc]
#Change format from string to DateTiem
df.dtutc=DateTime.(df.dtutc, "yyyy-mm-dd HH:MM:SS")
#Add BE to sample dictionary
push!(scenario_data["Generation"]["RES"]["Solar PV"]["BE"],"2020"=>df[1:8784,[:dtutc,:BE_solar_fc]])
push!(scenario_data["Generation"]["RES"]["Solar PV"]["BE"],"2021"=>df[8785:end,[:dtutc,:BE_solar_fc]])
push!(scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"],"2020"=>df[1:8784,[:dtutc,:BE_onshore_wind_fc]])
push!(scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"],"2021"=>df[8785:end,[:dtutc,:BE_onshore_wind_fc]])
#change columns names for BE
rename!(scenario_data["Generation"]["RES"]["Solar PV"]["BE"]["2020"], Symbol.(["time_stamp","BE_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"]["2020"], Symbol.(["time_stamp","BE_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Solar PV"]["BE"]["2021"], Symbol.(["time_stamp","BE_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"]["2021"], Symbol.(["time_stamp","BE_MWh"]))

#Add UK to sample dictionary
push!(scenario_data["Generation"]["RES"]["Solar PV"]["UK"],"2020"=>df[1:8784,[:dtutc,:UK_solar_fc]])
push!(scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"],"2020"=>df[1:8784,[:dtutc,:UK_wind_onshore_fc]])
push!(scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"],"2020"=>df[1:8784,[:dtutc,:UK_wind_offshore_fc]])
push!(scenario_data["Generation"]["RES"]["Solar PV"]["UK"],"2021"=>df[8785:end,[:dtutc,:UK_solar_fc]])
push!(scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"],"2021"=>df[8785:end,[:dtutc,:UK_wind_onshore_fc]])
push!(scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"],"2021"=>df[8785:end,[:dtutc,:UK_wind_offshore_fc]])
#change columns names for UK
rename!(scenario_data["Generation"]["RES"]["Solar PV"]["UK"]["2020"], Symbol.(["time_stamp","UK_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"]["2020"], Symbol.(["time_stamp","UK_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"]["2020"], Symbol.(["time_stamp","UK_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Solar PV"]["UK"]["2021"], Symbol.(["time_stamp","UK_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"]["2021"], Symbol.(["time_stamp","UK_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"]["2021"], Symbol.(["time_stamp","UK_MWh"]))

#do same for demand dataframe
scenario_data["Demand"]["Base"]["2020"]=DataFrame(time_stamp=df.dtutc[1:8784], BE_MWh=Float64.(df.BE_load_fc[1:8784]), UK_MWh=Float64.(df.UK_load_fc[1:8784]))
push!(scenario_data["Demand"]["NT"],"2021"=>DataFrame(time_stamp=df.dtutc[8785:end], BE_MWh=Float64.(df.BE_load_fc[8785:end]), UK_MWh=Float64.(df.UK_load_fc[8785:end])))
#load CorRES data
#df2=CSV.read("C://Users//shardy//Documents//julia//times_series_input_large_files//YUSO_data//P_DA_FC_35.csv", DataFrame)
df2=CSV.read("C://Users//shardy//Documents//julia//times_series_input_large_files//YUSO_data//P_DA_FC_58.csv", DataFrame)
rename!(df2,[:time,:BE_MWh])
df2[!,:time]=DateTime.(df2[!,:time], "mm/dd/yyyy HH:MM")
push!(scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"],"2020"=>df2[1:8784,[:time,:BE_MWh]])
push!(scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"],"2021"=>df2[8785:end,[:time,:BE_MWh]])
rename!(scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"]["2020"], Symbol.(["time_stamp","BE_MWh"]))
rename!(scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"]["2021"], Symbol.(["time_stamp","BE_MWh"]))
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBE58_2020.jld2",scenario_data)


#swap 2020 for 2021
scenario_data["Generation"]["RES"]["Solar PV"]["BE"]["2021"][!,:time_stamp]=scenario_data["Generation"]["RES"]["Solar PV"]["BE"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Generation"]["RES"]["Solar PV"]["BE"]["2020"]=scenario_data["Generation"]["RES"]["Solar PV"]["BE"]["2021"]

scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"]["2021"][!,:time_stamp]=scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"]["2020"]=scenario_data["Generation"]["RES"]["Onshore Wind"]["BE"]["2021"]

scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"]["2021"][!,:time_stamp]=scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"]["2020"]=scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"]["2021"]

scenario_data["Generation"]["RES"]["Solar PV"]["UK"]["2021"][!,:time_stamp]=scenario_data["Generation"]["RES"]["Solar PV"]["UK"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Generation"]["RES"]["Solar PV"]["UK"]["2020"]=scenario_data["Generation"]["RES"]["Solar PV"]["UK"]["2021"]

scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"]["2021"][!,:time_stamp]=scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"]["2020"]=scenario_data["Generation"]["RES"]["Onshore Wind"]["UK"]["2021"]

scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"]["2021"][!,:time_stamp]=scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"]["2020"]=scenario_data["Generation"]["RES"]["Offshore Wind"]["UK"]["2021"]

scenario_data["Demand"]["NT"]["2021"][!,:time_stamp]=scenario_data["Demand"]["NT"]["2021"][!,:time_stamp].-Year(1)
scenario_data["Demand"]["Base"]["2020"]=scenario_data["Demand"]["NT"]["2021"]
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBE58_2021.jld2",scenario_data)


################################################
#Processing from Hakan's Input file directly
################################################
scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKFRBENLDEDKNO.jld2")
scenario_data["Generation"]["Scenarios"]["Base"]["2020"]
scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4UKFRBENLDEDKNO.jld2")

scenario_data_4=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4UKFRBENLDEDKNO.jld2")


FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKFRBENLDEDKNO.jld2",scenario_data)


scenario_data["Generation"]["Scenarios"]=Dict()
push!(scenario_data["Generation"]["Scenarios"],"Base"=>Dict())
#push!(scenario_data["Generation"]["Scenarios"]["Base"],"2020"=>sdgs["NT2025"])

#Demand
sdgs=deepcopy(scenario_data_4["Demand"])
scenario_data["Demand"]=Dict()
push!(scenario_data["Demand"],"NT"=>Dict())
push!(scenario_data["Demand"],"GA"=>Dict())
push!(scenario_data["Demand"],"DE"=>Dict())
push!(scenario_data["Demand"],"Base"=>Dict())

push!(scenario_data["Demand"]["Base"],"2020"=>sdgs["NT2025"])
push!(scenario_data["Demand"]["NT"],"2030"=>sdgs["NT2030"])
push!(scenario_data["Demand"]["GA"],"2030"=>sdgs["GA2030"])
push!(scenario_data["Demand"]["DE"],"2030"=>sdgs["DE2030"])
push!(scenario_data["Demand"]["NT"],"2040"=>sdgs["NT2040"])
push!(scenario_data["Demand"]["GA"],"2040"=>sdgs["GA2040"])
push!(scenario_data["Demand"]["DE"],"2040"=>sdgs["DE2040"])

#scenario_data["Demand"]==scenario_data_for["Demand"]


#Generation
sdgs=deepcopy(scenario_data_4["Generation"])
scenario_data["Generation"]=Dict()
push!(scenario_data["Generation"],"keys"=>scenario_data_4["Generation"]["keys"])
push!(scenario_data["Generation"],"RES"=>scenario_data_4["Generation"]["RES"])
push!(scenario_data["Generation"],"costs"=>scenario_data_4["Generation"]["costs"])

push!(scenario_data["Generation"],"Scenarios"=>Dict())
push!(scenario_data["Generation"]["Scenarios"],"NT"=>Dict())
push!(scenario_data["Generation"]["Scenarios"],"GA"=>Dict())
push!(scenario_data["Generation"]["Scenarios"],"DE"=>Dict())
push!(scenario_data["Generation"]["Scenarios"],"Base"=>Dict())

push!(scenario_data["Generation"]["Scenarios"]["Base"],"2020"=>sdgs["Scenarios"]["NT2025"])
push!(scenario_data["Generation"]["Scenarios"]["NT"],"2030"=>sdgs["Scenarios"]["NT2030"])
push!(scenario_data["Generation"]["Scenarios"]["GA"],"2030"=>sdgs["Scenarios"]["GA2030"])
push!(scenario_data["Generation"]["Scenarios"]["DE"],"2030"=>sdgs["Scenarios"]["DE2030"])
push!(scenario_data["Generation"]["Scenarios"]["NT"],"2040"=>sdgs["Scenarios"]["NT2040"])
push!(scenario_data["Generation"]["Scenarios"]["GA"],"2040"=>sdgs["Scenarios"]["GA2040"])
push!(scenario_data["Generation"]["Scenarios"]["DE"],"2040"=>sdgs["Scenarios"]["DE2040"])



scenario_data["Generation"]["Scenarios"]["Base"]["2020"]["BLNK"]