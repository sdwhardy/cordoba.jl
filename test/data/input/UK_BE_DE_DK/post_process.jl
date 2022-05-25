################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

results_nodal=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\nodal_market_k6.jld2")
results_allhm=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_allhm_k6.jld2")
results_567=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_567_k6.jld2")

s_nodal=results_nodal["s"];result_mip_nodal=results_nodal["result_mip"];data_nodal=results_nodal["data"];mn_data_nodal=results_nodal["mn_data"]
s_allhm=results_allhm["s"];result_mip_allhm=results_allhm["result_mip"];data_allhm=results_allhm["data"];mn_data_allhm=results_allhm["mn_data"]
s_567=results_567["s"];result_mip_567=results_567["result_mip"];data_567=results_567["data"];mn_data_567=results_567["mn_data"]

s_nodal=_CBD.owpps_profit_obz(s_nodal, result_mip_nodal, mn_data_nodal)
s_allhm=_CBD.owpps_profit_obz(s_allhm, result_mip_allhm, mn_data_allhm)
s_567=_CBD.owpps_profit_obz(s_567, result_mip_567, mn_data_567)

s_nodal=_CBD.transmission_lines_profits(s_nodal, result_mip_nodal, mn_data_nodal, data_nodal);
s_allhm=_CBD.transmission_lines_profits(s_allhm, result_mip_allhm, mn_data_allhm, data_allhm);
s_567=_CBD.transmission_lines_profits(s_567, result_mip_567, mn_data_567, data_567);

s_nodal=_CBD.undo_marginal_price_scaling(s_nodal,result_mip_nodal)
s_allhm=_CBD.undo_marginal_price_scaling(s_allhm,result_mip_allhm)
s_567=_CBD.undo_marginal_price_scaling(s_567,result_mip_567)

gen_consume_summary_nodal=_CBD.summarize_generator_solution_data(result_mip_nodal, data_nodal,s_nodal)#print solution
gen_consume_summary_allhm=_CBD.summarize_generator_solution_data(result_mip_allhm, data_allhm,s_allhm)#print solution
gen_consume_summary_567=_CBD.summarize_generator_solution_data(result_mip_567, data_567,s_567)#print solution

social_welfare_nodal = _CBD.SocialWelfare(s_nodal, result_mip_nodal, mn_data_nodal, data_nodal)
social_welfare_allhm = _CBD.SocialWelfare(s_allhm, result_mip_allhm, mn_data_allhm, data_allhm)
social_welfare_567 = _CBD.SocialWelfare(s_567, result_mip_567, mn_data_567, data_567)

s["income_summary"]["owpp"]
social_welfare["totals"]
_CBD.topology_map(s,"tinf")
#nodal "gross_consumer_surplus"=>-1.30748e6
#hm
_CBD.plot_cumulative_production_all_scenarios_allWF(s, mn_data)
_CBD.plot_cumulative_income_all_scenarios_allWF(s, mn_data)
_CBD.plot_cumulative_income_tl_all_scenarios(s,data)

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
