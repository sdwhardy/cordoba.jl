################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

results_nodal=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\nodal_market_k6_VOLL5000b.jld2")
results_allhm=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_allhm_k6_VOLL5000b_rc.jld2")
results_567=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_567_k6_VOLL5000b_rc.jld2")




s_nodal, result_mip_nodal, data_nodal, mn_data_nodal=_CBD.summarize_in_s(results_nodal);
s_allhm, result_mip_allhm, data_allhm, mn_data_allhm=_CBD.summarize_zonal_in_s(results_allhm);
s_567, result_mip_567, data_567, mn_data_567=_CBD.summarize_zonal_in_s(results_567);

_CBD.print_table_summary(s_nodal)
_CBD.print_table_summary(s_allhm)
_CBD.print_table_summary(s_567)


_CBD.topology_map(s_567)

_CBD.plot_cumulative_production_all_scenarios_allWF(s_nodal, mn_data_nodal)
_CBD.plot_cumulative_income_all_scenarios_allWF(s_nodal, mn_data_nodal)
_CBD.plot_cumulative_income_tl_all_scenarios(s_nodal,data_nodal)


_CBD.print_solution_wcost_data(result_mip_allhm, s_allhm, data_allhm)

for (k,c) in s["income_summary"]["tso"]["totals"]["ac"];println("AC("*k*"):"*string(c));end

for (k,c) in s_allhm["income_summary"]["tso"]["totals"]["dc"];println("DC("*k*"):"*string(c));end

for (n,e) in result_mip_allhm["solution"]["nw"]
    dif=e["bus"]["1"]["lam_kcl_r"]*-1-e["bus"]["5"]["lam_kcl_r"]*-1
    if (e["branchdc"]["7"]["pf"]<0 && dif!=0)
    println(string(e["branchdc"]["7"]["pf"])*" "*string(dif))
end;end

country="UK"#,"DE","DK"]
scenario="3"
con=s_allhm["gen_consume_summary"]["onshore_demand"][scenario][country]#[121:144,:]
gen=s_allhm["gen_consume_summary"]["onshore_generation"][scenario][country]#[121:144,:]
_CBD.plot_generation_profile(deepcopy(gen),deepcopy(con),country*" "*scenario)

######################################################################
s_nodal=results_nodal["s"];result_mip_nodal=results_nodal["result_mip"];data_nodal=results_nodal["data"];mn_data_nodal=results_nodal["mn_data"]
s_allhm=results_allhm["s"];result_mip_allhm=results_allhm["result_mip"];data_allhm=results_allhm["data"];mn_data_allhm=results_allhm["mn_data"]
s_567=results_567["s"];result_mip_567=results_567["result_mip"];data_567=results_567["data"];mn_data_567=results_567["mn_data"]

s_nodal=_CBD.owpps_profit_obz(s_nodal, result_mip_nodal, mn_data_nodal)
s_allhm=_CBD.owpps_profit_obz(s_allhm, result_mip_allhm, mn_data_allhm)
s_567=_CBD.owpps_profit_obz(s_567, result_mip_567, mn_data_567)

s_nodal=_CBD.transmission_lines_profits(s_nodal, result_mip_nodal, mn_data_nodal, data_nodal);
s_allhm=_CBD.transmission_lines_profits(s_allhm, result_mip_allhm, mn_data_allhm, data_allhm);
s_567=_CBD.transmission_lines_profits(s_567, result_mip_567, mn_data_567, data_567);

s_nodal=_CBD.undo_marginal_price_scaling(s_nodal, result_mip_nodal)
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
_CBD.plot_cumulative_income_all_scenarios_allWF(s_567, mn_data_567)
_CBD.plot_cumulative_income_tl_all_scenarios(s_allhm,data_allhm)

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

country="UK"#,"DE","DK"]
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

##################### nOBZ in zonal config
###########################################
#NOTE uncomment AC cables!!!!!!!!!!!!!!!!! 1043 post_process.jl
s=results_nodal["s"]
result_mip=results_nodal["result_mip"]

mn_data=results_nodal["mn_data"]

data=results_nodal["data"]
zones=[[2,5],[3,6],[4,7]]
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
s["home_market"]=zones
s["rebalancing"]=true
s["relax_problem"]=true
s["output"]["duals"]=true
result_mip_hm_prices = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_wcost_data(result_mip_hm_prices, s, data)#-856896.0245340846 

s["home_market"]=[]
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
result_mip= _CBD.hm_market_prices(result_mip, result_mip_hm_prices)
_CBD.print_solution_wcost_data(result_mip, s, data)#-856896.0245340846 
results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s, "result_mip_hm_prices"=>result_mip_hm_prices)

s, result_mip, data, mn_data=_CBD.summarize_zonal_in_s(results);

results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s)

s, result_mip, data, mn_data=_CBD.summarize_in_s(results);

_CBD.print_table_summary(s)
###########################################
#################### printong 2 at once #############################
node="5";scenario="1"
time_series_wfuk=s_nodal["income_summary"]["onshore"][scenario]["1"]["price_not_npv"]
time_series_wfbe=s_nodal["income_summary"]["onshore"][scenario]["2"]["price_not_npv"]
time_series_wfde=s_nodal["income_summary"]["onshore"][scenario]["3"]["price_not_npv"]
time_series_wfdk=s_nodal["income_summary"]["onshore"][scenario]["4"]["price_not_npv"]
for scenario in ["2","3","4","5","6"]
time_series_wfbe=time_series_wfbe.+s_nodal["income_summary"]["onshore"][scenario]["2"]["price_not_npv"]
time_series_wfde=time_series_wfde.+s_nodal["income_summary"]["onshore"][scenario]["3"]["price_not_npv"]
time_series_wfdk=time_series_wfdk.+s_nodal["income_summary"]["onshore"][scenario]["4"]["price_not_npv"] 
end
time_series_wfbe=time_series_wfbe/7
time_series_wfde=time_series_wfde/7
time_series_wfdk=time_series_wfdk/7

function plot_clearing_price(time_series)
    
    clrs=generation_color_map()
    

        #low_rng=minimum(marginal_prices)
        #high_rng=maximum(marginal_prices)
        scatter_vec_gen=[
            
            PlotlyJS.scatter(
                y=time_series_wfbe,
                mode="lines", name="BE",
                line=PlotlyJS.attr(width=2, color="red")
            ),PlotlyJS.scatter(
                y=time_series_wfde,
                mode="lines", name="DE",
                line=PlotlyJS.attr(width=2, color="blue")
            ),PlotlyJS.scatter(
                y=time_series_wfdk,
                mode="lines", name="DK",
                line=PlotlyJS.attr(width=2, color="green")
            ) ]
        PlotlyJS.plot(
            scatter_vec_gen, PlotlyJS.Layout(font_size=35,xaxis_range=(0, 432),yaxis_title="€/MWh",xaxis_title="time steps"))
    end
