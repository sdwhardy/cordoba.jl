################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

#################################### map output #######################
############################## newest #################################
results_NT=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\NORTH_SEA_nodal_k4_NT.jld2")
results_DE=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\NORTH_SEA_nodal_k4_DE.jld2")
results_GA=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\NORTH_SEA_nodal_k4_GA.jld2")

result=results_GA
pdic=_CBD.problemOUTPUT_map_byTimeStep(result["6"])
PlotlyJS.plot(pdic["trace012"], pdic["layout"])

for (i, res) in result
    for (j,dc_br) in res["result_mip"]["solution"]["nw"]["576"]["branchdc"]
        if (dc_br["p_rateA"]>0.0)
            f_bus=0;t_bus=0;
            for (k,conv) in result[i]["data"]["convdc"]
                if (conv["busdc_i"] == res["data"]["branchdc"][j]["fbusdc"]) 
                    f_bus=conv["busac_i"];end
                if (conv["busdc_i"] == res["data"]["branchdc"][j]["tbusdc"]) 
                    t_bus=conv["busac_i"];end
            end
            println(j,": from ",f_bus," to ",t_bus," MVA: ",dc_br["p_rateA"])
        end
    end
end

for (i, res) in result
    for (j,dc_br) in res["result_mip"]["solution"]["nw"]["576"]["branch"]
        if (dc_br["p_rateAC"]>0.0)
        
            f_bus=res["data"]["branch"][j]["f_bus"]
                
            t_bus=res["data"]["branch"][j]["t_bus"]

            println(j,": from ",f_bus," to ",t_bus," MVA: ",dc_br["p_rateAC"])
        end
    end
end

######################################################################################################
results=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\NORTH_SEA_nodal_k4_full.jld2")

pdic=_CBD.problemOUTPUT_map_byTimeStep(result["1"])
PlotlyJS.plot(pdic["trace012"], pdic["layout"])

resultMIP=deepcopy(result)

resultMIP["1"]["result_mip"]=result["1"]["result_mip2"]["solution"]["nw"]
pdicmip=_CBD.problemMIP_OUTPUT_map_byTimeStep(resultMIP["1"])
PlotlyJS.plot(pdicmip["trace012"], pdicmip["layout"])

result["1"]["data"]["branchdc"]["42"]["branchdc"]
for i=1:1:1728 
println(i," ", result["1"]["result_mip"]["solution"]["nw"][string(i)]["branchdc"]["42"]["p_rateA"]," ", result["1"]["result_mip"]["solution"]["nw"][string(i)]["branchdc"]["42"]["pt"])
end
i=1728
for i=1:1:1728 
    println(i," ", result["1"]["result_mip"]["solution"]["nw"][string(i)]["gen"]["14"]["pg"])
end

for i=1:1:1728 
    println(i," ", result["1"]["result_mip"]["solution"]["nw"][string(i)]["convdc"]["29"]["p_pacmax"])
end
i=1
cbdc=[]
for i=1:1:1728
    for (c_i, cdc) in result["1"]["result_mip"]["solution"]["nw"][string(i)]["branchdc"]
        if (result["1"]["data"]["branchdc"][c_i]["fbusdc"]>31 && result["1"]["data"]["branchdc"][c_i]["tbusdc"]>31)
            if (isapprox(cdc["p_rateA"], cdc["pt"], atol=0.01))
                #println(c_i, " ", cdc["p_rateA"], " ", cdc["pt"], " ", result["1"]["data"]["branchdc"][c_i]["fbusdc"], " ", result["1"]["data"]["branchdc"][c_i]["tbusdc"])
                push!(cbdc, c_i)
            end
        end
    end
end


for c_i in cbdc
    fac=0;tac=0
    for (cv_i,cv) in result["1"]["data"]["convdc"]
        if (result["1"]["data"]["branchdc"][c_i]["fbusdc"]==cv["busdc_i"])
            fac=cv["busac_i"]
        end
        if (result["1"]["data"]["branchdc"][c_i]["tbusdc"]==cv["busdc_i"])
            tac=cv["busac_i"]
        end
    end
    println(c_i, " ", fac, " ", tac)
end


unique!(cbdc)