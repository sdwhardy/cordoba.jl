################## loads external packages ##############################
using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames, Clp
using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
import cordoba; const _CBD = cordoba#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
include("../aux/post_process/functions.jl")

############################# data cleaning script
################### ENTSO-E scenario description ####################################
scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
scenario_years=["2020","2030","2040"]

#k=100
#scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series

#scenario_data = Dict{String,Any}()

byr="17"#17,18,19,20
ryr="2020"#2020,2030,2040
_sc="EU";#EU,ST,DG
for byr in ["20"]#["17","18","19","20"]#
for ryr in ["2020","2030","2040"]
for _sc in ["EU","ST","DG"]
        if (ryr=="2020")
                n="0"
        elseif (_sc=="EU" && ryr=="2030")
                n="1"
        elseif (_sc=="ST" && ryr=="2030")
                n="2"
        elseif (_sc=="DG" && ryr=="2030")
                n="3"
        elseif (_sc=="EU" && ryr=="2040")
                n="4"
        elseif (_sc=="ST" && ryr=="2040")
                n="5"
        elseif (_sc=="DG" && ryr=="2040")
                n="6"
        end

        #loads wnd series and handles exceptions
        df_be=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/BE_"*n*".csv", DataFrames.DataFrame)
        rename!(df_be, Dict(Symbol("Day-ahead Price [EUR/MWh]") => "EUR_daBE")); if (haskey(df_be,"Column1"));select!(df_be, Not(:Column1));end
        if (haskey(df_be,"BZN|BE"));select!(df_be, Not(Symbol("BZN|BE")));end
        df_de=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/DE_"*n*".csv", DataFrames.DataFrame)
        rename!(df_de, Dict(Symbol("Day-ahead Price [EUR/MWh]") => "EUR_daDE")); if (haskey(df_de,"Column1"));select!(df_de, Not(:Column1));end
        df_dk=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/DK_"*n*".csv", DataFrames.DataFrame)
        rename!(df_dk, Dict(Symbol("Day-ahead Price [EUR/MWh]") => "EUR_daDK")); if (haskey(df_dk,"Column1"));select!(df_dk, Not(:Column1));end
        df_uk=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/UK_"*n*".csv", DataFrames.DataFrame)
        rename!(df_uk, Dict(Symbol("Day-ahead Price [EUR/MWh]") => "EUR_daUK")); if (haskey(df_uk,"Column1"));select!(df_uk, Not(:Column1));end
        if (haskey(df_uk,"Day-ahead Price [GBP/MWh]"));select!(df_uk, Not(Symbol("Day-ahead Price [GBP/MWh]")));end
        if (haskey(df_uk,"Day-ahead Price [BGP/MWh]"));select!(df_uk, Not(Symbol("Day-ahead Price [BGP/MWh]")));end

        if (typeof(df_uk["EUR_daUK"][1])==String);floatz=[];for (i,ps) in  enumerate(df_uk["EUR_daUK"])
                        if (ps=="#VALUE!");push!(floatz,missing);else;push!(floatz,ps)
                        end;end;df_uk=hcat(df_uk,floatz);
                        df_uk=df_uk[completecases(df_uk), :]
                        if (haskey(df_uk,"EUR_daUK"));select!(df_uk, Not(Symbol("EUR_daUK")));end;
                        rename!(df_uk, Dict(Symbol("x1") => "EUR_daUK"));
                        df_uk["EUR_daUK"]=[parse(Float64,i) for i in df_uk["EUR_daUK"]]
                end


        df_price=innerjoin(df_be,df_de,on = Symbol("MTU (CET)"), makeunique=true)
        df_price=innerjoin(df_price,df_dk,on = Symbol("MTU (CET)"), makeunique=true)
        df_price=innerjoin(df_price,df_uk,on = Symbol("MTU (CET)"), makeunique=true)
        dt_price=[];for i=1:1:length(df_price["MTU (CET)"]);push!(dt_price,DateTime(df_price["MTU (CET)"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
        df_price=hcat(df_price,dt_price); if (haskey(df_price,"MTU (CET)"));select!(df_price, Not(Symbol("MTU (CET)")));end;
        rename!(df_price, Dict(Symbol("x1") => "time_stamp"))

        df_wBE=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/BE_wnd.csv", DataFrames.DataFrame)
        dt_be=[];for i=1:1:length(df_wBE["MTU"]);push!(dt_be,DateTime(df_wBE["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
        df_wBE=hcat(df_wBE,dt_be); if (haskey(df_wBE,"MTU"));select!(df_wBE, Not(Symbol("MTU")));end; if (haskey(df_wBE,"Area"));select!(df_wBE, Not(Symbol("Area")));end
        rename!(df_wBE, Dict(Symbol("x1") => "time_stamp"));rename!(df_wBE, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]") => "Wnd_MWhBE"));

        df_wDE=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/DE_wnd.csv", DataFrames.DataFrame)
        dt_de=[];for i=1:1:length(df_wDE["MTU"]);push!(dt_de,DateTime(df_wDE["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
        df_wDE=hcat(df_wDE,dt_de); if (haskey(df_wDE,"MTU"));select!(df_wDE, Not(Symbol("MTU")));end; if (haskey(df_wDE,"Area"));select!(df_wDE, Not(Symbol("Area")));end
        rename!(df_wDE, Dict(Symbol("x1") => "time_stamp"));rename!(df_wDE, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]") => "Wnd_MWhDE"));

        df_wDK=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/DK_wnd.csv", DataFrames.DataFrame)
        dt_dk=[];for i=1:1:length(df_wDK["MTU"]);push!(dt_dk,DateTime(df_wDK["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
        df_wDK=hcat(df_wDK,dt_dk); if (haskey(df_wDK,"MTU"));select!(df_wDK, Not(Symbol("MTU")));end; if (haskey(df_wDK,"Area"));select!(df_wDK, Not(Symbol("Area")));end
        rename!(df_wDK, Dict(Symbol("x1") => "time_stamp"));rename!(df_wDK, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]") => "Wnd_MWhDK"));

        df_wUK=CSV.read("./test/data/input/scenarios/time series data_w_BE/Based_on_20"*byr*"_prices/UK_wnd.csv", DataFrames.DataFrame)
        dt_uk=[];for i=1:1:length(df_wUK["MTU"]);push!(dt_uk,DateTime(df_wUK["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
        df_wUK=hcat(df_wUK,dt_uk); if (haskey(df_wUK,"MTU"));select!(df_wUK, Not(Symbol("MTU")));end; if (haskey(df_wUK,"Area"));select!(df_wUK, Not(Symbol("Area")));end
        rename!(df_wUK, Dict(Symbol("x1") => "time_stamp"));rename!(df_wUK, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]") => "Wnd_MWhUK"));

        df_wnd=innerjoin(df_wBE,df_wDE,on = Symbol("time_stamp"), makeunique=true)
        df_wnd=innerjoin(df_wnd,df_wDK,on = Symbol("time_stamp"), makeunique=true)
        df_wnd=innerjoin(df_wnd,df_wUK,on = Symbol("time_stamp"), makeunique=true)

        df_wnd=innerjoin(df_wnd,df_price,on = Symbol("time_stamp"), makeunique=true)

        df_wnd=df_wnd[completecases(df_wnd), :]
        df_wnd=filter(:Wnd_MWhUK => x -> x!="N/A", df_wnd)
        df_wnd=filter(:Wnd_MWhUK => x -> !any(f -> f(x), (ismissing, isnothing)), df_wnd)
        #df_wnd["Wnd_MWhUK"]=[parse(Int64,i) for i in df_wnd["Wnd_MWhUK"]]
        df_wnd["Wnd_MWhBE"]=df_wnd["Wnd_MWhBE"]./maximum(df_wnd["Wnd_MWhBE"])
        df_wnd["Wnd_MWhDK"]=df_wnd["Wnd_MWhDK"]./maximum(df_wnd["Wnd_MWhDK"])
        df_wnd["Wnd_MWhUK"]=df_wnd["Wnd_MWhUK"]./maximum(df_wnd["Wnd_MWhUK"])
        if (haskey(df_wnd,"Currency"));select!(df_wnd, Not(Symbol("Currency")));end;
        if (haskey(df_wnd,"Currency_1"));select!(df_wnd, Not(Symbol("Currency_1")));end;
        if (haskey(df_wnd,"Currency_2"));select!(df_wnd, Not(Symbol("Currency_2")));end;
        display(df_wnd)
        CSV.write("./test/data/input/scenarios/convex_wBE/"*_sc*byr*ryr*".csv",df_wnd)
end;end;end

###########################################################################
##################### Cluster time series data ###########################
################### ENTSO-E scenario description ####################################
scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
#scenario_names=["EU17","EU18","EU19","EU20"]
#scenario_names=["EU19"]

scenario_years=["2020","2030","2040"]
#scenario_years=["2020","2030"]
#scenario_years=["2030"]

##################### load time series data ##############################
k=100
scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series
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
        daily_bew=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhBE"]))
        daily_dkw=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhDK"]))
        daily_ukw=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhUK"]))
        daily_be=_CBD.daily_tss(Float64.(ts[!,"EUR_daBE"]))
        daily_de=_CBD.daily_tss(Float64.(ts[!,"EUR_daDE"]))
        daily_uk=_CBD.daily_tss(Float64.(ts[!,"EUR_daUK"]))
        daily_dk=_CBD.daily_tss(Float64.(ts[!,"EUR_daDK"]))

        kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_ukw.^2).+(daily_dkw.^2).+(daily_bew.^2).+(daily_dew.^2).+(daily_be.^2).+(daily_de.^2).+(daily_dk.^2).+(daily_uk.^2))',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
        filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
    end
end
save("./test/data/input/BE_DE_UK_DK_islands/EUSTDG_17to20_TS_k"*string(k)*"_wBE.jld2",scenario_data)


########################################## Depricated
#=if (haskey(df,"MTU (UTC)"));rename!(df, Dict(Symbol("MTU (UTC)") => "MTU"));end
dt=[];for i=1:1:length(df["MTU"]);push!(dt,DateTime(df["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
dt_w=[];for i=1:1:length(df_w["MTU"]);push!(dt_w,DateTime(df_w["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
dt_wDK=[];for i=1:1:length(df_wDK["MTU"]);push!(dt_wDK,DateTime(df_wDK["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
dt_wUK=[];for i=1:1:length(df_wUK["MTU"]);push!(dt_wUK,DateTime(df_wUK["MTU"][i][1:16],dateformat"dd.mm.yyyy HH:MM"));end
df=hcat(df,dt);df_w=hcat(df_w,dt_w);df_wDK=hcat(df_wDK,dt_wDK);df_wUK=hcat(df_wUK,dt_wUK);df_temp=hcat(df,df_w,makeunique=true);
df_wUK=df_wUK[completecases(df_wUK), :]
df_temp=hcat(df_temp,df_wDK,makeunique=true);df=DataFrame()
rename!(df_temp, Dict(Symbol("x1") => "time_stamp"))
rename!(df_wUK, Dict(Symbol("x1") => "time_stamp"))
df_temp=innerjoin(df_temp,df_wUK,on = :time_stamp, makeunique=true)
df=innerjoin(df_temp,orig,on = :time_stamp)
if (haskey(df,"x1_2"));select!(df, Not(:x1_2));end
if (haskey(df,"x1"));select!(df, Not(:x1));end
if (haskey(df,"x1_1"));select!(df, Not(:x1_1));end
if (haskey(df,"Wnd_MWhDK"));select!(df, Not(:Wnd_MWhDK));end
if (haskey(df,"Wnd_MWhUK"));select!(df, Not(:Wnd_MWhUK));end
if (haskey(df,"BZN|BE"));select!(df, Not(Symbol("BZN|BE")));end
if (haskey(df,"Column1"));select!(df, Not(Symbol("Column1")));end
if (haskey(df,"MTU (UTC)"));select!(df, Not(Symbol("MTU (UTC)")));end
if (haskey(df,"MTU"));select!(df, Not(Symbol("MTU")));end
if (haskey(df,"MTU_1"));select!(df, Not(Symbol("MTU_1")));end
if (haskey(df,"MTU_2"));select!(df, Not(Symbol("MTU_2")));end
if (haskey(df,"Area"));select!(df, Not(Symbol("Area")));end
if (haskey(df,"Area_1"));select!(df, Not(Symbol("Area_1")));end
rename!(df, Dict(Symbol("Day-ahead Price [EUR/MWh]") => "EUR_daBE"))
rename!(df, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]") => "Wnd_MWhBE"))
rename!(df, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]_1") => "Wnd_MWhDK"))
rename!(df, Dict(Symbol("Wind Offshore  - Actual Aggregated [MW]_2") => "Wnd_MWhUK"))



byr="19"#17,18,19,20
ryr="2020"#2020,2030,2040
_sc="EU";#EU,ST,DG
for byr in ["20"]#["17","18","19","20"]
for ryr in ["2020","2030","2040"]
for _sc in ["EU","ST","DG"]
        if (ryr=="2020")
                n="0"
        elseif (_sc=="EU" && ryr=="2030")
                n="1"
        elseif (_sc=="ST" && ryr=="2030")
                n="2"
        elseif (_sc=="DG" && ryr=="2030")
                n="3"
        elseif (_sc=="EU" && ryr=="2040")
                n="4"
        elseif (_sc=="ST" && ryr=="2040")
                n="5"
        elseif (_sc=="DG" && ryr=="2040")
                n="6"
        end
        df=CSV.read("./test/data/input/scenarios/convex_wBE/"*_sc*byr*ryr*".csv", DataFrames.DataFrame)

        floatz=[]
        if (typeof(df["EUR_daUK"][1])==String)
                for (i,ps) in  enumerate(scenario_data["EU19"]["2020"]["EUR_daUK"])
                        println(i)
                        if (ps=="#VALUE!")
                                push!(floatz,missing)
                        else
                                push!(floatz,parse(Float64,ps))
                        end
                end
        end
        scenario_data["EU19"]["2020"]["EUR_daUK"]=#
