
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
    end
    return extradata
end


function add_storage_profile(dim, data, extradata, zs_data, zs, number_of_hours)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU

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
        #up_reg
        up_pe=[(zs_data["EUR_up"*z][d],zs_data["MWh_up"*z][d]) for z in zs]
        sort!(up_pe, by = x -> x[1], rev=true)
        up_pe=(first.(up_pe)./e2me,last.(up_pe)./pu)

        dwn_pe=[(zs_data["EUR_dwn"*z][d],zs_data["MWh_dwn"*z][d]) for z in zs]
        elec_prices=[(zs_data["EUR_id"*z][d],Inf) for z in zs]
        sort!(elec_prices, by = x -> x[1])
        #push!(dwn_pe,elec_prices[1])
        sort!(dwn_pe, by = x -> x[1])
        dwn_pe=(first.(dwn_pe)./e2me,last.(dwn_pe)./pu)


        for (b,bat) in storage
            #Set discharge rating and cost
            dcr_temp=deepcopy(bat["discharge_rating"])
            extradata["ne_storage"][b]["discharge_rating"][1, d] = 0
            extradata["ne_storage"][b]["cost_inj"][1, d] = 0
            for upe=1:1:length(up_pe[1])
                if (dcr_temp<=last(up_pe)[upe])
                    extradata["ne_storage"][b]["discharge_rating"][1, d] = extradata["ne_storage"][b]["discharge_rating"][1, d] + dcr_temp
                    extradata["ne_storage"][b]["cost_inj"][1, d] = extradata["ne_storage"][b]["cost_inj"][1, d] + first(up_pe)[upe]
                    last(up_pe)[upe]=last(up_pe)[upe]-dcr_temp
                    break
                elseif (last(up_pe)[upe]>1e-3 && last(up_pe)[upe]<dcr_temp)
                    extradata["ne_storage"][b]["discharge_rating"][1, d] = extradata["ne_storage"][b]["discharge_rating"][1, d] + last(up_pe)[upe]
                    extradata["ne_storage"][b]["cost_inj"][1, d] = ((extradata["ne_storage"][b]["discharge_rating"][1, d]-last(up_pe)[upe])/extradata["ne_storage"][b]["discharge_rating"][1, d])*extradata["ne_storage"][b]["cost_inj"][1, d] + (last(up_pe)[upe]/extradata["ne_storage"][b]["discharge_rating"][1, d])*first(up_pe)[upe]
                    dcr_temp=dcr_temp-last(up_pe)[upe]
                    last(up_pe)[upe] = 0
                elseif (last(up_pe)[upe]<=1e-3 && extradata["ne_storage"][b]["discharge_rating"][1, d] <= 1e-3)
                    extradata["ne_storage"][b]["discharge_rating"][1, d] = deepcopy(bat["discharge_rating"])
                    extradata["ne_storage"][b]["cost_inj"][1, d] = 0
                    break
                end
            end
            #Set charge rating and cost
            cr_temp=deepcopy(bat["charge_rating"])
            extradata["ne_storage"][b]["charge_rating"][1, d] = 0
            extradata["ne_storage"][b]["cost_abs"][1, d] = 0
            for upe=1:1:length(dwn_pe[1])
                if (cr_temp<=last(dwn_pe)[upe] && last(dwn_pe)[upe]!=Inf)
                    extradata["ne_storage"][b]["charge_rating"][1, d] = extradata["ne_storage"][b]["charge_rating"][1, d] + cr_temp
                    extradata["ne_storage"][b]["cost_abs"][1, d] = extradata["ne_storage"][b]["cost_abs"][1, d] + first(dwn_pe)[upe]
                    last(dwn_pe)[upe]=last(dwn_pe)[upe]-cr_temp
                    break
                elseif (last(dwn_pe)[upe]>1e-3 && last(dwn_pe)[upe]<cr_temp)
                    extradata["ne_storage"][b]["charge_rating"][1, d] = extradata["ne_storage"][b]["charge_rating"][1, d] + last(dwn_pe)[upe]
                    extradata["ne_storage"][b]["cost_abs"][1, d] = ((extradata["ne_storage"][b]["charge_rating"][1, d]-last(dwn_pe)[upe])/extradata["ne_storage"][b]["charge_rating"][1, d])*extradata["ne_storage"][b]["cost_abs"][1, d] + (last(dwn_pe)[upe]/extradata["ne_storage"][b]["charge_rating"][1, d])*first(dwn_pe)[upe]
                    cr_temp=cr_temp-last(dwn_pe)[upe]
                    last(dwn_pe)[upe] = 0
                elseif ((last(dwn_pe)[upe]<=1e-3 || last(dwn_pe)[upe]==Inf) && extradata["ne_storage"][b]["charge_rating"][1, d] <= 1e-3)
                    extradata["ne_storage"][b]["charge_rating"][1, d] = deepcopy(bat["charge_rating"])
                    extradata["ne_storage"][b]["cost_abs"][1, d] = first(dwn_pe)[upe]
                    break
                end
            end
        end
    end

    return extradata,data
end

function create_profile_sets_wstrg(number_of_hours, data, df0, df1,ic_mva,owpp_mva)
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

function create_profile_sets_mesh(number_of_hours, data_orig, zs_data, zs, inf_grid, owpp_mva)
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    data=Dict{String,Any}();data["gen"]=Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    data["gen"]=sort(data_orig["gen"])
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end
    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
        for (g, gen) in data["gen"]
            if (gen["apf"]>0)#market generator onshore
                extradata["gen"][g]["pmax"][1, d] = inf_grid/pu
                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [(zs_data["EUR_da"*zs[gen["gen_bus"]]][d])/e2me,0]
            else#wind gen
                extradata["gen"][g]["pmax"][1, d] = (zs_data["Wnd_MWh"*zs[gen["gen_bus"]]][d])*owpp_mva/pu
                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [0,0]
            end
        end
    end
    #add loads
    loads=Dict{String,Any}()
    num_of_gens=length(data["gen"])
    for (g, gen) in extradata["gen"]
        if (data["gen"][g]["apf"]>0)#market generator onshore
            load=deepcopy(data["gen"][g])
            load["index"]=num_of_gens+1
            load["source_id"][2]=num_of_gens+1
            load["pmin"]=deepcopy(load["pmax"])*-1
            load["pmax"]=0
            push!(loads,string(num_of_gens+1)=>deepcopy(load))
            num_of_gens=num_of_gens+1
        end
    end
    for (l, load) in loads
        extradata["gen"][l] = Dict{String,Any}()
        extradata["gen"][l]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][l]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][l]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
        for (l, load) in loads
            if (load["apf"]>0)#market generator onshore
                extradata["gen"][l]["pmax"][1, d] = 0
                extradata["gen"][l]["pmin"][1, d] = (inf_grid/pu)*-1
                extradata["gen"][l]["cost"][d] = [(zs_data["EUR_da"*zs[load["gen_bus"]]][d])/e2me,0]
                push!(data_orig["gen"],l=>load)
            else#wind gen
            end
        end
    end

    #set ["apf"]
    for (g, gen) in data_orig["gen"]
        gen["apf"]=0
    end
    return extradata,data_orig
end

function get_profile_data_sets_mesh(zs,data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    zs_data=[]
    for z in zs
        df=CSV.read("./test/data/input/"*z*"data.csv", DataFrames.DataFrame)
        colnames = ["time_stamp","Wnd_MWh"*z,"EUR_da"*z,"EUR_id"*z,"MWh_up"*z,"EUR_up"*z,"MWh_dwn"*z,"EUR_dwn"*z]
        names!(df, Symbol.(colnames))
        push!(zs_data,df)
    end
    zsd=zs_data[1];for z in zs_data[2:end];
    zsd=innerjoin(zsd,z, makeunique=true,on=:time_stamp);end
    Ytr=PCA_cluster(zsd)
    Ytr = convert(Matrix, zsd[1:1:end,2:end])'
    cluster=kmeans_cluster(Ytr)
    #cluster=kmedoid_cluster(Ytr)
    n_samples=n_samps(cluster,n)
    zsd=zsd[n_samples,:];
    s="1";scenario["sc_years"][s]
    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]
        scenario["hours"]=length(zsd[!,:time_stamp])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, zsd
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
