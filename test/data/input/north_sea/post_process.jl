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