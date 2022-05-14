################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

results=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\nodal_test_run.jld2")
s=results["s"];result_mip=results["result_mip"];data=results["data"];mn_data=results["mn_data"]
#_CBD.print_solution_wcost_data(result_mip, s, data)#-856896.0245340846
s=_CBD.owpps_profit_obz(s, result_mip, mn_data)
s=_CBD.transmission_lines_profits(s, result_mip, mn_data, data);


_CBD.topology_map(s,"tinf")
   
_CBD.plot_cumulative_wf_income_all_scenarios(s, mn_data, "DE")

_CBD.plot_cumulative_wf_income_all_scenarios(s, mn_data, "DK")

_CBD.plot_cumulative_wf_income_all_scenarios(s, mn_data, "BE")

_CBD.plot_cumulative_income_tl_all_scenarios(s,data)




c=sum(sum(cb["rent"]) for (c,cb) in hourly_income_tl["dc"])#3262.67155
c=sum(sum(cb["rent"]) for (c,cb) in hourly_income_tl["ac"])



country="BE"#,"DE","DK"]
scenario="1"
con=gen_consume_summary["onshore_demand"][scenario][country]#[121:144,:]
gen=gen_consume_summary["onshore_generation"][scenario][country]#[121:144,:]


#result_mip=deepcopy(result_mip_001)
gen_consume_summary=_CBD.summarize_generator_solution_data(result_mip, data,s)#print solution
#####################################

country="DE"#,"DE","DK"]
scenario="1"
con=gen_consume_summary["onshore_demand"][scenario][country]#[121:144,:]
gen=gen_consume_summary["onshore_generation"][scenario][country]#[121:144,:]
_CBD.plot_generation_profile(deepcopy(gen),deepcopy(con),country*" "*scenario)


s["income_summary"]=Dict()

country="DE"#,"DE","DK"]
scenario="1"
con=gen_consume_summary["onshore_demand"][scenario][country]#[121:144,:]
gen=gen_consume_summary["offshore_generation"][scenario][country]#[121:144,:]
con=select!(con,:ts)
_CBD.plot_generation_profile(deepcopy(gen),deepcopy(con),country*" "*scenario)


result_mip=_CBD.undo_marginal_price_scaling(s,result_mip)
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (1,"UK"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (2,"BE"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (3,"DE"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (4,"DK"))

social_welfare = SocialWelfare(s, result_mip, mn_data, data)
social_welfare["totals"]
