function create_profile_sets(number_of_hours, data, df0, df1,ic_mva,owpp_mva)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU
    #e2me=1
    da=0.835;id=0.165
    #e2me=1#into ME/PU
    extradata = Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    storage=[(i,b) for (i,b) in data["ne_storage"]]
    sort!(storage, by=x->x[2]["energy_rating"])
    extradata["ne_storage"] = Dict{String,Any}()
    for (b, bat) in data["ne_storage"]
        extradata["ne_storage"][b] = Dict{String,Any}()
        extradata["ne_storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["charge_rating"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["discharge_rating"] = Array{Float64,2}(undef, 1, number_of_hours)
    end

    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df0.wind[d]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [da*df0.daprice[d]/e2me+id*df0.idprice[d]/e2me,0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df0.wind[d]
            extradata["gen"]["2"]["cost"][d] = [da*df0.daprice[d]/e2me+id*df0.idprice[d]/e2me,0]

            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df0.wind[d]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [da*df1.daprice[d]/e2me+id*df1.idprice[d]/e2me,0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df0.wind[d]
                extradata["gen"]["4"]["cost"][d] = [da*df1.daprice[d]/e2me+id*df1.idprice[d]/e2me,0]
            #Wind generator
                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/pu*df0.wind[d]
                extradata["gen"]["5"]["pmin"][1, d] = 0
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]

                #up_reg
                up_prices=[df0.regup_price[d],df1.regup_price[d]];
                best_up_price=deepcopy(findmax(up_prices));up_prices[best_up_price[2]]=Inf
                worst_up_price=findmin(up_prices)
                best_up_mwh=[df0.regup_mwh[d],df1.regup_mwh[d]][best_up_price[2]]
                worst_up_mwh=[df0.regup_mwh[d],df1.regup_mwh[d]][worst_up_price[2]]
                #down_reg
                dwn_prices=[df0.regdwn_price[d],df1.regdwn_price[d]];
                best_dwn_price=deepcopy(findmax(dwn_prices));dwn_prices[best_dwn_price[2]]=Inf
                worst_dwn_price=findmin(dwn_prices)
                best_dwn_mwh=[df0.regdwn_mwh[d],df1.regdwn_mwh[d]][best_dwn_price[2]]
                worst_dwn_mwh=[df0.regdwn_mwh[d],df1.regdwn_mwh[d]][worst_dwn_price[2]]

                best_up_mwh = best_up_mwh>0 ? best_up_mwh/pu : 0
                best_dwn_mwh = best_dwn_mwh>0 ? best_dwn_mwh/pu : 0
                worst_up_mwh = worst_up_mwh>0 ? worst_up_mwh/pu : 0
                worst_dwn_mwh = worst_dwn_mwh>0 ? worst_dwn_mwh/pu : 0

                best_up_price = best_up_price[1]>0 ? best_up_price[1]/e2me : 0
                best_dwn_price = best_dwn_price[1]>0 ? best_dwn_price[1]/e2me : 0
                worst_up_price = worst_up_price[1]>0 ? worst_up_price[1]/e2me : 0
                worst_dwn_price = worst_dwn_price[1]>0 ? worst_dwn_price[1]/e2me : 0

                for (b,bat) in storage
                    dcr_temp=deepcopy(data["ne_storage"][b]["discharge_rating"])
                    bup_mwh_temp=deepcopy(best_up_mwh)
                    wdwn_mwh_temp=deepcopy(worst_up_mwh)

                    best_up_mwh = dcr_temp>best_up_mwh ?  best_up_mwh : dcr_temp#can it be from the battery?
                    best_up_mwh = best_up_mwh>0 ?  best_up_mwh : 0
                    worst_up_mwh = dcr_temp>(best_up_mwh+worst_up_mwh) ?  worst_up_mwh : dcr_temp-best_up_mwh#can additional be from the battery?
                    worst_up_mwh = worst_up_mwh>0 ?  worst_up_mwh : 0
                    if ((best_up_mwh+worst_up_mwh)>0.0001)
                        extradata["ne_storage"][b]["discharge_rating"][1, d]=best_up_mwh+worst_up_mwh
                        extradata["ne_storage"][b]["cost_inj"][1, d] = (best_up_mwh/(best_up_mwh+worst_up_mwh))*best_up_price+((worst_up_mwh)/(best_up_mwh+worst_up_mwh))*worst_up_price
                        best_up_mwh=bup_mwh_temp-best_up_mwh
                        worst_up_mwh=wdwn_mwh_temp-worst_up_mwh
                    else
                        extradata["ne_storage"][b]["discharge_rating"][1, d]=0
                        extradata["ne_storage"][b]["cost_inj"][1, d]=0
                    end
                    if (isnan(extradata["ne_storage"][b]["discharge_rating"][1, d]))
                        println("Bad data - detected at "*string(d)*" - "*string(b))
                        extradata["ne_storage"][b]["discharge_rating"][1, d]=0
                        extradata["ne_storage"][b]["cost_inj"][1, d]=0
                    end
                    cr_temp=deepcopy(data["ne_storage"][b]["charge_rating"])
                    bdwn_mwh_temp=deepcopy(best_dwn_mwh)
                    wdwn_mwh_temp=deepcopy(worst_dwn_mwh)

                    best_dwn_mwh = cr_temp>best_dwn_mwh ?  best_dwn_mwh : cr_temp#can it be from the battery?
                    best_dwn_mwh = best_dwn_mwh>0 ?  best_dwn_mwh : 0
                    worst_dwn_mwh = cr_temp>(best_dwn_mwh+worst_dwn_mwh) ?  worst_dwn_mwh : cr_temp-best_dwn_mwh#can additional be from the battery?
                    worst_dwn_mwh = worst_dwn_mwh>0 ?  worst_dwn_mwh : 0
                    if ((best_dwn_mwh+worst_dwn_mwh)>0.0001)
                        extradata["ne_storage"][b]["charge_rating"][1, d]=cr_temp
                        extradata["ne_storage"][b]["cost_abs"][1, d] = (best_dwn_mwh/cr_temp)*best_dwn_price+(worst_dwn_mwh/cr_temp)*worst_dwn_price
                        best_dwn_mwh=bdwn_mwh_temp-best_dwn_mwh
                        worst_dwn_mwh=wdwn_mwh_temp-worst_dwn_mwh
                    else
                        extradata["ne_storage"][b]["charge_rating"][1, d]=cr_temp
                        extradata["ne_storage"][b]["cost_abs"][1, d]=0
                    end
                end
    end
    return extradata
end

function scale_bat_data_cordoba!(data, scenario)
    rescale_hourly = x -> (scenario["hours"] / (8760*scenario["planning_horizon"])) * x # yearly limit on energy absoption
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "max_energy_absorption", rescale_hourly)
    end
end

function scale_cost_data_cordoba!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / scenario["hours"]) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / scenario["hours"]) * x # scale total costs to the planning horizon
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
        _PM._apply_func!(branch, "co2_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
        _PM._apply_func!(branch, "co2_cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
        _PM._apply_func!(conv, "co2_cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "co2_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end
end


function get_profile_data_sets(d1,d2,data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];

    df0=DataFrame("time_stamp"=>[],"daprice"=>[],"idprice"=>[],"wind"=>[],"regup_price"=>[],"regup_mwh"=>[],"regdwn_price"=>[],"regdwn_mwh"=>[])
    df1=DataFrame("time_stamp"=>[],"daprice"=>[],"idprice"=>[],"wind"=>[],"regup_price"=>[],"regup_mwh"=>[],"regdwn_price"=>[],"regdwn_mwh"=>[])#MWh_up,EUR_up,MWh_dwn,EUR_dwn
    z0_data=CSV.read("./test/data/input/"*d1*".csv", DataFrames.DataFrame)
    z1_data=CSV.read("./test/data/input/"*d2*".csv", DataFrames.DataFrame)
    z01_data=innerjoin(z0_data,z1_data, makeunique=true,on=:time_stamp)
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]

        #tss=[DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.be_eumwh.+0.165.*ukbe_data.be_costid),DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.uk_eumwh.+0.165.*ukbe_data.uk_costid)]
        tss=[DataFrame("time_stamp"=>z01_data.time_stamp, "price"=>abs.((0.835.*z01_data.EUR_da.+0.165.*z01_data.EUR_id).-(0.835.*z01_data.EUR_da_1.+0.165.*z01_data.EUR_id_1))),DataFrame("time_stamp"=>z01_data.time_stamp, "price"=>z01_data.Wnd_MWh)]
        tss_bins=cluster_ts(tss,n)
        sc=sample_cluster(tss_bins,tss,n)
        sort!(sc)
        for t in sc
            #Zone 0
            push!(df0,[t,z0_data[!,:EUR_da][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:EUR_id][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:Wnd_MWh][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:MWh_up][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:EUR_up][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:MWh_dwn][findfirst(isequal(t),z0_data[!,:time_stamp])],
            z0_data[!,:EUR_dwn][findfirst(isequal(t),z0_data[!,:time_stamp])]])
            #Zone 1
            push!(df1,[t,z1_data[!,:EUR_da][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:EUR_id][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:Wnd_MWh][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:MWh_up][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:EUR_up][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:MWh_dwn][findfirst(isequal(t),z1_data[!,:time_stamp])],
            z1_data[!,:EUR_dwn][findfirst(isequal(t),z1_data[!,:time_stamp])]])
        end

        scenario["hours"]=length(df0[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, df0,df1
end


function get_n_profile_data(data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    if haskey(scenario, "mc")
        monte_carlo = scenario["mc"]
    else
        monte_carlo = false
    end
    df=DataFrame("time_stamp"=>[],"z0_daprice"=>[],"z0_idprice"=>[],"z1_daprice"=>[],"z1_idprice"=>[],"z0_wind"=>[],"z1_wind"=>[])
    ukbe_data=CSV.read("./test/data/cordoba/input/ukbe_ts_2.csv", DataFrames.DataFrame)
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]

        #tss=[DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.be_eumwh.+0.165.*ukbe_data.be_costid),DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>0.835.*ukbe_data.uk_eumwh.+0.165.*ukbe_data.uk_costid)]
        tss=[DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>abs.((0.835.*ukbe_data.be_eumwh.+0.165.*ukbe_data.be_costid).-(0.835.*ukbe_data.uk_eumwh.+0.165.*ukbe_data.uk_costid))),DataFrame("time_stamp"=>ukbe_data.time_stamp, "price"=>ukbe_data.be_wind)]
        tss_bins=cluster_ts(tss,n)
        sc=sample_cluster(tss_bins,tss,n)
        sort!(sc)
        for t in sc
            push!(df,[t,ukbe_data[!,:be_eumwh][findfirst(isequal(t),ukbe_data[!,:time_stamp])], ukbe_data[!,:be_costid][findfirst(isequal(t),ukbe_data[!,:time_stamp])],ukbe_data[!,:uk_eumwh][findfirst(isequal(t),ukbe_data[!,:time_stamp])], ukbe_data[!,:uk_costid][findfirst(isequal(t),ukbe_data[!,:time_stamp])],ukbe_data[!,:be_wind][findfirst(isequal(t),ukbe_data[!,:time_stamp])], ukbe_data[!,:uk_wind][findfirst(isequal(t),ukbe_data[!,:time_stamp])]])
        end

        scenario["hours"]=length(df[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, df
end
