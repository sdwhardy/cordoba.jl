#"connection,cost,time,ic_mva,ic_km,owpp_mva,owpp_km"
function store_data(wcfile,resultACDC, ic_mva, owpp_mva, ic_length, owpp_km)
    if (resultACDC["solution"]["nw"]["1"]["ne_branch"]["1"]["built"]>0)
        connection=-1
    elseif (resultACDC["solution"]["nw"]["1"]["ne_branch"]["2"]["built"]>0)
        connection=0
    else
        connection=1
    end
    println(wcfile,string(connection)*", "*string(resultACDC["objective"])*", "*string(resultACDC["solve_time"])*", "*string(ic_mva)*", "*string(ic_length)*", "*string(owpp_mva)*", "*string(owpp_km))
end


# Example code to print list of built HVDC branches and converters
function display_results(result)
    built_cv = []
    built_br = []
    built_ACbr = []
    for (c, conv) in result["solution"]["convdc_ne"]
        if isapprox(conv["isbuilt"] , 1; atol = 0.01)
            print("Conv: $c \n")
            push!(built_cv,c)
        end
    end
    for (b, branch) in result["solution"]["branchdc_ne"]
        if isapprox(branch["isbuilt"] , 1; atol = 0.01)
            print("DCBranch: $b \n")
            push!(built_br,b)
        end
    end
    for (b, branch) in result["solution"]["ne_branch"]
        if isapprox(branch["built"] , 1; atol = 0.01)
            print("ACBranch: $b \n")
            push!(built_ACbr,b)
        end
    end
    return built_cv, built_br, built_ACbr
end
