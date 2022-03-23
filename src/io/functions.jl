
function rebuild_system_converters(konverters,data_struct, off_shore_nodes,cost)
    for (key,cnv) in konverters
        data_struct["convdc"][key]["Pacmin"]=-1*cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Pacmax"]=cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Qacmin"]=-1*cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Qacmax"]=cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Pacrated"]=cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Qacrated"]=cnv["p_pacmax"]*100
        if (issubset(parse(Int64,key),off_shore_nodes))
            data_struct["convdc"][key]["cost"]=0.29*cnv["p_pacmax"]*100
        else
            data_struct["convdc"][key]["cost"]=0.1925*cnv["p_pacmax"]*100
        end
        cost=cost+data_struct["convdc"][key]["cost"]
    end
    return data_struct,cost
end

function rebuild_system_cables_convex(kable_indices,data,cost,z_base)
    for k in kable_indices
        if (first(last(k))>0)
            cb=_CBD.DC_cbl(first(last(k)),last(last(k)))
            data["branchdc"][first(k)]["rateA"]=data["branchdc"][first(k)]["rateB"]=data["branchdc"][first(k)]["rateC"]=cb.num*cb.elec.mva
            data["branchdc"][first(k)]["cost"]=cb.costs.cpx_i+cb.costs.cpx_p
            data["branchdc"][first(k)]["r"]=((cb.elec.ohm*10^3/cb.num)*cb.length)/z_base
            cost=cost+data["branchdc"][first(k)]["cost"]
        else
            data["branchdc"][first(k)]["rateA"]=data["branchdc"][first(k)]["rateB"]=data["branchdc"][first(k)]["rateC"]=0
            data["branchdc"][first(k)]["cost"]=0
        end
    end
    return data,cost
end

function rebuild_system_cables(kables,data,cost,hpy)
    for (key,brn) in kables
        if (haskey(brn, "rateA") && brn["rateA"]>0)
            data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=brn["rateA"]*100
            data["branchdc"][key]["cost"]=brn["cost"]*hpy
            data["branchdc"][key]["r"]=brn["r"]
            cost=cost+data["branchdc"][key]["cost"]
        elseif (haskey(brn, "rateA"))
            data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=0
            data["branchdc"][key]["cost"]=0
        elseif (haskey(brn, "rateA") && brn["rateA"]>0)
            data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=brn["rateA"]*100
            data["branchdc"][key]["cost"]=brn["cost"]*hpy
            data["branchdc"][key]["r"]=brn["r"]
            cost=cost+data["branchdc"][key]["cost"]
        elseif (haskey(brn, "rateA"))
            data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=0
            data["branchdc"][key]["cost"]=0
        end;
    end
    return data,cost
end


function print_topology_data_DC(data_mip,markets_wfs)
    println("%%%%%%%%%%%%%%%%%%%%%%%%% Nodes %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,n) in enumerate(markets_wfs[1])
        println(string(i)*" - "*n)
    end
    for (i,n) in enumerate(markets_wfs[2])
        println(string(length(markets_wfs[1])+i)*" - OWPP_"*n)
    end
    println("%%%%%%%%%%%%%%%%%% Converters %%%%%%%%%%%%%%%%%%%%")
    for (i,cv) in sort(OrderedCollections.OrderedDict(data_mip["convdc"]), by=x->parse(Int64,x))
            println(string(i)*": "*string(cv["cost"]))
    end;
    println("%%%%%%%%%%%%%%%%%%%%%%%% Cables %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,br) in sort(OrderedCollections.OrderedDict(data_mip["branchdc_ne"]), by=x->parse(Int64,x))
        println(string(i)*": "*string(br["fbusdc"])*" - "*string(br["tbusdc"])*" MVA: "*string(br["rateA"])*" Length: "*string(br["length"])*" cost: "*string(br["cost"]))
    end
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
end

function print_topology_data_AC(data_mip,markets_wfs)
    println("%%%%%%%%%%%%%%%%%%%%%%%%% Nodes %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,n) in enumerate(markets_wfs[1])
        println(string(i)*" - "*n)
    end
    for (i,n) in enumerate(markets_wfs[2])
        println(string(length(markets_wfs[1])+i)*" - OWPP_"*n)
    end
    println("%%%%%%%%%%%%%%%%%%%%%%%% Cables %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,br) in sort(OrderedCollections.OrderedDict(data_mip["ne_branch"]), by=x->parse(Int64,x))
        println(string(i)*": "*string(br["f_bus"])*" - "*string(br["t_bus"])*" MVA: "*string(br["rate_a"])*" Length: "*string(br["length"])*" Cost: "*string(br["construction_cost"]))
    end
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
end



function print_solution_data(result_mip, data_mip, argz)
    println("Description: test-"*string(argz["test"])*" k-"*string(argz["k"])*" years-"*string(argz["scenario_years"])*" scenarios-"*string(argz["scenario_names"]))
    if (haskey(result_mip["solution"]["nw"]["1"],"branch"))
        println("%%%%%%% CONVEX SOLUTION %%%%%%%")
        print_branch(result_mip,argz,data_mip)
        print_branchdc(result_mip,argz,data_mip)
    else
        println("%%%%%%% MIP SOLUTION %%%%%%%")
        print_branch_ne(result_mip,argz,data_mip)
        print_branchdc_ne(result_mip,argz,data_mip)
    end
    print_owpps(result_mip,argz)
    print_converters(result_mip,argz)
    print_storage(result_mip,argz)
    println("objective: "*string(result_mip["objective"])*" achieved in: "*string(result_mip["solve_time"]))
end

#=function print_convex_solution_data(result_mip, data_mip, argz)
    print_owpps(result_mip,argz)
    print_branch(result_mip,argz,data_mip)
    print_branchdc(result_mip,argz,data_mip)
    print_converters(result_mip,argz)
    print_storage(result_mip,argz)
    println("objective: "*string(result_mip["objective"]))
end=#

function print_branch(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"branch"))
        println("%%%% Cables HVAC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["branch"]), by=x->parse(Int64,x))
            println(string(i)*": "*string(data_mip["branch"][i]["f_bus"])*" - "*string(data_mip["branch"][i]["t_bus"])*" MVA: "*string(br["p_rateAC"]))
        end
        println("%%%% Cables HVAC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["branch"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branch"][i]["f_bus"])*" - "*string(data_mip["branch"][i]["t_bus"])*" MVA: "*string(br["p_rateAC"]))
        end
        println("%%%% Cables HVAC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branch"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branch"][i]["f_bus"])*" - "*string(data_mip["branch"][i]["t_bus"])*" MVA: "*string(br["p_rateAC"]))
        end
    end
end

function print_branchdc(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"branchdc"))
        println("%%%% Cables HVDC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["branchdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branchdc"][i]["fbusdc"])*" - "*string(data_mip["branchdc"][i]["tbusdc"])*" MVA: "*string(br["p_rateA"]))
        end
        println("%%%% Cables HVDC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["branchdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branchdc"][i]["fbusdc"])*" - "*string(data_mip["branchdc"][i]["tbusdc"])*" MVA: "*string(br["p_rateA"]))
        end
        println("%%%% Cables HVDC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branchdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(data_mip["branchdc"][i]["fbusdc"])*" - "*string(data_mip["branchdc"][i]["tbusdc"])*" MVA: "*string(br["p_rateA"]))
        end
    end
end

function print_branch_ne(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"ne_branch"))
        println("%%%% Cables HVAC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["ne_branch"]), by=x->parse(Int64,x))
            if (br["built"]==1)
                println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
        end;end
        println("%%%% Cables HVAC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["ne_branch"]), by=x->parse(Int64,x))
            if (br["built"]==1)
                println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
        end;end
        println("%%%% Cables HVAC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["ne_branch"]), by=x->parse(Int64,x))
            if (br["built"]==1)
                println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
        end;end
    end
end

function print_branchdc_ne(result_mip,argz,data_mip)
    if (haskey(result_mip["solution"]["nw"]["1"],"branchdc_ne"))
        println("%%%% Cables HVDC t0 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["branchdc_ne"]), by=x->parse(Int64,x))
            if (br["isbuilt"]==1)
                println(string(i)*": "*string(data_mip["branchdc_ne"][i]["fbusdc"])*" - "*string(data_mip["branchdc_ne"][i]["tbusdc"])*" MVA: "*string(data_mip["branchdc_ne"][i]["rateA"])*" cost: "*string(data_mip["branchdc_ne"][i]["cost"]))
        end;end
        println("%%%% Cables HVDC t2 %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["branchdc_ne"]), by=x->parse(Int64,x))
            if (br["isbuilt"]==1)
                println(string(i)*": "*string(data_mip["branchdc_ne"][i]["fbusdc"])*" - "*string(data_mip["branchdc_ne"][i]["tbusdc"])*" MVA: "*string(data_mip["branchdc_ne"][i]["rateA"])*" cost: "*string(data_mip["branchdc_ne"][i]["cost"]))
        end;end
        println("%%%% Cables HVDC tinf %%%%")
        for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["branchdc_ne"]), by=x->parse(Int64,x))
            if (br["isbuilt"]==1)
                println(string(i)*": "*string(data_mip["branchdc_ne"][i]["fbusdc"])*" - "*string(data_mip["branchdc_ne"][i]["tbusdc"])*" MVA: "*string(data_mip["branchdc_ne"][i]["rateA"])*" cost: "*string(data_mip["branchdc_ne"][i]["cost"]))
        end;end
    end
end

function print_converters(result_mip,argz)
    if (haskey(result_mip["solution"]["nw"]["1"],"convdc"))
        println("%%%% Converters t0 %%%%")
        for (i,cv) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["convdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(cv["p_pacmax"]))
        end;
        println("%%%% Converters t2 %%%%")
        for (i,cv) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["convdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(cv["p_pacmax"]))
        end;
        println("%%%% Converters tinf %%%%")
        for (i,cv) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["convdc"]), by=x->parse(Int64,x))
                println(string(i)*": "*string(cv["p_pacmax"]))
        end;
    end
end

function print_storage(result_mip,argz)
    if (haskey(result_mip["solution"]["nw"]["1"],"storage"))
        println("%%%% Storage t0 %%%%")
        for (i,s) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["storage"]), by=x->parse(Int64,x))
                println(string(i)*": "*" MWh: "*string(s["e_absmax"]))
        end
        println("%%%% Storage t2 %%%%")
        for (i,s) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["storage"]), by=x->parse(Int64,x))
                println(string(i)*": "*" MWh: "*string(s["e_absmax"]))
        end
        println("%%%% Storage tinf %%%%")
        for (i,s) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["storage"]), by=x->parse(Int64,x))
                println(string(i)*": "*" MWh: "*string(s["e_absmax"]))
        end
    end
end

function print_owpps(result_mip,argz)
    if (haskey(result_mip["solution"]["nw"]["1"],"gen"))
        println("%%%% OWPPS T0 %%%%")
        for (i,wf) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]["1"]["gen"]), by=x->parse(Int64,x))
            if (haskey(wf,"wf_pacmax"))
                println(string(i)*": "*string(wf["wf_pacmax"]))
        end;end
        println("%%%% OWPPS T2 %%%%")
        for (i,wf) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(argz["ls"]+1)]["gen"]), by=x->parse(Int64,x))
            if (haskey(wf,"wf_pacmax"))
                println(string(i)*": "*string(wf["wf_pacmax"]))
        end;end
        println("%%%% OWPPS Tinf %%%%")
        for (i,wf) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["gen"]), by=x->parse(Int64,x))
            if (haskey(wf,"wf_pacmax"))
                println(string(i)*": "*string(wf["wf_pacmax"]))
        end;end
    end
end


function print_solution_data_AC(result_mip, data_mip)
    println("%%%%%%%%%%%%%%%%%%%%%%%% Cables %%%%%%%%%%%%%%%%%%%%%%%")
    for (i,br) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"][string(length(result_mip["solution"]["nw"]))]["ne_branch"]), by=x->parse(Int64,x))
        if (br["built"]==1)
            println(string(i)*": "*string(data_mip["ne_branch"][i]["f_bus"])*" - "*string(data_mip["ne_branch"][i]["t_bus"])*" MVA: "*string(data_mip["ne_branch"][i]["rate_a"])*" cost: "*string(data_mip["ne_branch"][i]["construction_cost"]))
    end;end

end
