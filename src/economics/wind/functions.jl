
################################################################################
######################### Wind Calculations ####################################
################################################################################
#**
function wndF_wndPrf(profiles)
    wp=Array{Float32,1}()
    dummy=Array{Float32,1}()
    dummy=profiles[1]
    for i=2:length(profiles)
        dummy=dummy+profiles[i]
    end
    wp=dummy./length(profiles)

    wnd=wind()
    mx=findmax(wp)[1]
    ord=reverse(sort(wp))./mx
    wndF_conEng(ord,wnd)#constraint energy calc
    wndF_ldLss(ord, wnd)#calculates loss factor
    wnd
    return wnd
end
#loss factor
#**
function wndF_ldLss(div, wind)
  #wind.lf=(sum(div))*0.85/length(div)#saves loss factor, 0.85 is wake penalization
  wind.lf=(sum(div))/length(div)#saves loss factor, 0.85 is wake penalization
  #loss factor/llf formula ref: Guidelines on the calculation and use of loss factors Te Mana Hiko Electricity Authority
  #0.85 for wake effect from Evaluation of the wind direction uncertainty and its impact on wake modeling at the Horns Rev offshore wind farm
  #M. Gaumond  P.‐E. Réthoré  S. Ott  A. Peña  A. Bechmann  K. S. Hansen
  llf=0.0
  for pu in div
    #llf=llf+(pu*0.85)^2
    llf=llf+pu^2
  end
  wind.delta=llf/length(div)#saves load loss factor
end

#Finds Constraind energy
#**
function  wndF_conEng(graph,wnd)
    #create sized arrays
    ce=Float32[]
    for hr=1:length(graph)
        smPu=0
        for pu=1:hr
            smPu=smPu+graph[pu]
        end
        smPu=smPu-(graph[hr]*hr)
        push!(ce,smPu)
    end
    wnd.ce=deepcopy(ce)
    wnd.pu=deepcopy(graph)
    return nothing
end

#check
#**Totals wind profiles of set
function find_netWind(wndstring)
    if !(haskey(wind_module.wind_profs,wndstring))
        wind_sum=wind()
        wndstring=wndstring[1]=='_' ? wndstring[2:end] : wndstring;
        wnds=split(wndstring,"_")
        if (length(wnds)>0 && wnds[length(wnds)]!="")
            wind_sum.pu=zeros(Float32,8759)
            wind_sum.ce=zeros(Float32,8759)
            wind_sum.delta=0
            wind_sum.lf=0
            true_length=0
            for w in wnds
                if (length(wind_module.wind_profs[w].pu)>0)
                    wind_sum.pu=(wind_sum.pu.+wind_module.wind_profs[w].pu)
                    wind_sum.ce=(wind_sum.ce.+wind_module.wind_profs[w].ce)
                    wind_sum.delta=(wind_sum.delta+wind_module.wind_profs[w].delta)
                    wind_sum.lf=(wind_sum.lf+wind_module.wind_profs[w].lf)
                    true_length=true_length+1
                end
            end
            if true_length!=0
                wind_sum.pu=wind_sum.pu./true_length
                wind_sum.ce=wind_sum.ce./true_length
                wind_sum.delta=wind_sum.delta/true_length
                wind_sum.lf=wind_sum.lf/true_length
            end
        end
        #println("storing: "*wndstring)
        wind_module.save_wind4_module(wind_sum,wndstring)
    else
    end
    return wndstring
end
