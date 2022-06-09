################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

results_nodal=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_DE_DK\\nodal_results_VOLL.jld2")
results_14=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_DE_DK\\zonal_results_hm14_VOLL.jld2")
results_24=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_DE_DK\\zonal_results_hm24_VOLL.jld2")
results_34=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_DE_DK\\zonal_results_hm34_VOLL5000.jld2")

s_nodal, result_mip_nodal, data_nodal, mn_data_nodal=_CBD.summarize_in_s(results_nodal);
s_14, result_mip_14, data_14, mn_data_14=_CBD.summarize_in_s(results_14);
s_24, result_mip_24, data_24, mn_data_24=_CBD.summarize_in_s(results_24);
s_34, result_mip_34, data_34, mn_data_34=_CBD.summarize_in_s(results_34);

_CBD.print_solution_wcost_data(result_mip_14, s_14, data_14)#-856896.0245340846
_CBD.print_table_summary(s_nodal)
_CBD.print_table_summary(s_14)
_CBD.print_table_summary(s_24)
_CBD.print_table_summary(s_34)


_CBD.topology_map(s_nodal,"tinf")

_CBD.plot_cumulative_production_all_scenarios_allWF(s_nodal, mn_data_nodal)
_CBD.plot_cumulative_income_all_scenarios_allWF(s_nodal, mn_data_nodal)
_CBD.plot_cumulative_income_tl_all_scenarios(s_nodal,data_nodal)












####################################################################################################################

s=s_nodal; result_mip=result_mip_nodal;mn_data=mn_data_nodal;data=data_nodal
sw=SocialWelfare(s_nodal, result_mip_nodal, mn_data_nodal, data_nodal)
s_nodal["totals"]["all"]
2.39119e8/2.84666e6
2.75453e8/2.84666e6
2.75453e8-2.39119e8

s_nodal=results_nodal["s"];result_mip_nodal=results_nodal["result_mip"];data_nodal=results_nodal["data"];mn_data_nodal=results_nodal["mn_data"]
s_14=results_14["s"];result_mip_14=results_14["result_mip"];data_14=results_14["data"];mn_data_14=results_14["mn_data"]
s_24=results_24["s"];result_mip_24=results_24["result_mip"];data_24=results_24["data"];mn_data_24=results_24["mn_data"]
s_34=results_34["s"];result_mip_34=results_34["result_mip"];data_34=results_34["data"];mn_data_34=results_34["mn_data"]

s_nodal=_CBD.owpps_profit_obz(s_nodal, result_mip_nodal, mn_data_nodal)
s_14=_CBD.owpps_profit_obz(s_14, result_mip_14, mn_data_14)
s_24=_CBD.owpps_profit_obz(s_24, result_mip_24, mn_data_24)
s_34=_CBD.owpps_profit_obz(s_34, result_mip_34, mn_data_34)

s_nodal=_CBD.transmission_lines_profits(s_nodal, result_mip_nodal, mn_data_nodal, data_nodal);
s_14=_CBD.transmission_lines_profits(s_14, result_mip_14, mn_data_14, data_14);
s_24=_CBD.transmission_lines_profits(s_24, result_mip_24, mn_data_24, data_24);
s_34=_CBD.transmission_lines_profits(s_34, result_mip_34, mn_data_34, data_34);

s_nodal=_CBD.undo_marginal_price_scaling(s_nodal,result_mip_nodal)
s_14=_CBD.undo_marginal_price_scaling(s_14,result_mip_14)
s_24=_CBD.undo_marginal_price_scaling(s_24,result_mip_24)
s_34=_CBD.undo_marginal_price_scaling(s_34,result_mip_34)

gen_consume_summary_nodal=_CBD.summarize_generator_solution_data(result_mip_nodal, data_nodal,s_nodal)#print solution
gen_consume_summary_14=_CBD.summarize_generator_solution_data(result_mip_14, data_14,s_14)#print solution
gen_consume_summary_24=_CBD.summarize_generator_solution_data(result_mip_24, data_24,s_24)#print solution
gen_consume_summary_34=_CBD.summarize_generator_solution_data(result_mip_34, data_34,s_34)#print solution


social_welfare_nodal = _CBD.SocialWelfare(s_nodal, result_mip_nodal, mn_data_nodal, data_nodal)
social_welfare_14 = _CBD.SocialWelfare(s_14, result_mip_14, mn_data_14, data_14)
social_welfare_24 = _CBD.SocialWelfare(s_24, result_mip_24, mn_data_24, data_24)
social_welfare_34 = _CBD.SocialWelfare(s_34, result_mip_34, mn_data_34, data_34)

s["income_summary"]["owpp"]
social_welfare["totals"]
_CBD.topology_map(s_nodal,"tinf")
#nodal "gross_consumer_surplus"=>-1.30748e6
#hm

s_34=_CBD.tl_totals(s_34,data_34)
s_34["income_summary"]["tso"]["totals"]
s_nodal["income_summary"]["strg"]["all"]


s=s_nodal;result_mip=result_mip_nodal;mn_data=mn_data_nodal
scenario_num="1";scenario= mn_data["scenario"][scenario_num]
s_nodal=_CBD.strg_profit_obzs(s_nodal, result_mip_nodal, mn_data_nodal)



_CBD.plot_cumulative_wf_production_all_scenarios(s, mn_data, "DE")

_CBD.plot_cumulative_wf_income_all_scenarios(s, mn_data, "DK")

_CBD.plot_cumulative_wf_income_all_scenarios(s, mn_data, "BE")

using PlotlyJS
yvalues=[social_welfare_allhm["totals"]["all"]["gross_consumer_surplus"], social_welfare_567["totals"]["all"]["gross_consumer_surplus"], social_welfare_nodal["totals"]["all"]["gross_consumer_surplus"]]
yvalues=[s_allhm["income_summary"]["owpp"]["all"]["power"]*(8760/s_allhm["hours_length"]), s_567["income_summary"]["owpp"]["all"]["power"]*(8760/s_567["hours_length"]), s_nodal["income_summary"]["owpp"]["all"]["power"]*(8760/s_nodal["hours_length"])]
xvalues=["HMD", "zOBZ", "nOBZ"]
xlabel="Market Type"
ylabel="GWh"
_CBD.basic_bar_chart(xvalues,yvalues,xlabel,ylabel)





#social welfare
#totals: -74468.4(567), -75017.9(HM), -74738.6(zonal)
c=sum(sum(cb["rent"]) for (c,cb) in hourly_income_tl["dc"])#3262.67155
c=sum(sum(cb["rent"]) for (c,cb) in hourly_income_tl["ac"])

wf=[];sc="6"
for (n,nw) in enumerate(s["xd"]["gen"]["218"]["pmax"]); if (issubset([n],values(mn_data["scenario"][sc]))) push!(wf,nw);end;end
maximum(wf)


country="BE"#,"DE","DK"]
scenario="1"
con=gen_consume_summary["onshore_demand"][scenario][country]#[121:144,:]
gen=gen_consume_summary["onshore_generation"][scenario][country]#[121:144,:]


#result_mip=deepcopy(result_mip_001)
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


_CBD.plot_dual_marginal_price(result_mip, string.(values(mn_data["scenario"][scenario])), (1,"UK"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (2,"BE"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (4,"BE"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (4,"DK"))
