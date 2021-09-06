
function storage_costs(data,year::Int64=2021)
    if (year==2021)
        #csts=[944,472,236,118,59,29.5,14.75,7.375]
        #csts=[472,236,118,59,29.5,14.75,7.375]
        csts=[236,118,59,29.5,14.75,7.375]
    elseif (year==2030)
        #csts=[560,280,140,70,35,17.5,8.75,4.375]
        #csts=[280,140,70,35,17.5,8.75,4.375]
        x=1
        #csts=[140/x,70/x,35/x,17.5/x,8.75/x,4.375/x]
        csts=[35/x,17.5/x,8.75/x,4.375/x,2.1875/x,1.09375/x]
        #csts=[1,1,1,1,1,1]
    else
        println("No battery cost data for specified year defaulting to 2021.")
        #csts=[944,472,236,118,59,29.5,14.75,7.375]
        #csts=[472,236,118,59,29.5,14.75,7.375]
        csts=[236,118,59,29.5,14.75,7.375]
    end
    [data["ne_storage"][string(i)]["inst_cost"]=data["ne_storage"][string(i)]["eq_cost"]=cst for (i,cst) in enumerate(csts)]
    [data["ne_storage"][string(i)]["eq_cost"]=cst for (i,cst) in enumerate(csts)]
end

function pm_dict_owpp(data,hoa)
    data["ne_branch"]["1"]=ne_branch(data["ne_branch"]["1"], hoa["candidate_lines"]["ac"]["13"][1], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_xfo(data["ne_branch"]["1"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_xfo(data["ne_branch"]["1"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_plat(data["ne_branch"]["1"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["ne_branch"]["2"]=ne_branch(data["ne_branch"]["2"], hoa["candidate_lines"]["ac"]["13"][2], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_xfo(data["ne_branch"]["2"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_xfo(data["ne_branch"]["2"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_plat(data["ne_branch"]["2"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["ne_branch"]["3"]=ne_branch(data["ne_branch"]["3"], hoa["candidate_lines"]["ac"]["13"][3], data["gen"]["1"]["mbase"])
    data["ne_branch"]["3"]=branch_xfo(data["ne_branch"]["3"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["3"]=branch_xfo(data["ne_branch"]["3"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["3"]=branch_plat(data["ne_branch"]["3"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["ne_branch"]["4"]=ne_branch(data["ne_branch"]["4"], hoa["candidate_lines"]["ac"]["13"][4], data["gen"]["1"]["mbase"])
    data["ne_branch"]["4"]=branch_xfo(data["ne_branch"]["4"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["4"]=branch_xfo(data["ne_branch"]["4"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["4"]=branch_plat(data["ne_branch"]["4"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["ne_branch"]["5"]=ne_branch(data["ne_branch"]["5"], hoa["candidate_lines"]["ac"]["13"][5], data["gen"]["1"]["mbase"])
    data["ne_branch"]["5"]=branch_xfo(data["ne_branch"]["5"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["5"]=branch_xfo(data["ne_branch"]["5"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["5"]=branch_plat(data["ne_branch"]["5"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["branchdc_ne"]["1"]=branchdc_ne(data["branchdc_ne"]["1"],hoa["candidate_lines"]["dc"]["12"],data["gen"]["1"]["mbase"])

    data["convdc_ne"]["1"]=convdc_ne(data["convdc_ne"]["1"],hoa["nodes"]["1"],data["gen"]["1"]["mbase"])
    data["convdc_ne"]["2"]=convdc_ne(data["convdc_ne"]["2"],hoa["nodes"]["2"],data["gen"]["1"]["mbase"])

    data["gen"]["1"]=gen(data["gen"]["1"],hoa["nodes"]["1"]["mva"]/data["gen"]["1"]["mbase"])
    data["gen"]["2"]=gen(data["gen"]["2"],-1*hoa["nodes"]["2"]["mva"]/data["gen"]["1"]["mbase"])
    data["gen"]["3"]=gen(data["gen"]["3"],hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["4"]=gen(data["gen"]["4"],-1*hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["5"]=gen(data["gen"]["5"],-1*hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    return data
end

function pm_dict_reducedBins(data,hoa)
    data["ne_branch"]["1"]=ne_branch(data["ne_branch"]["1"], hoa["candidate_lines"]["ac"]["13"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_xfo(data["ne_branch"]["1"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_xfo(data["ne_branch"]["1"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_plat(data["ne_branch"]["1"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["ne_branch"]["2"]=ne_branch(data["ne_branch"]["2"], hoa["candidate_lines"]["ac"]["23"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_xfo(data["ne_branch"]["2"],hoa["nodes"]["2"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_xfo(data["ne_branch"]["2"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_plat(data["ne_branch"]["2"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    #data["branchdc_ne"]["1"]=branchdc_ne(data["branchdc_ne"]["1"],hoa["candidate_lines"]["dc"]["12"],data["gen"]["1"]["mbase"])
    data["branchdc"]["1"]=branchdc(data["branchdc"]["1"],hoa["candidate_lines"]["dc"]["13"],data["gen"]["1"]["mbase"])
    data["branchdc"]["2"]=branchdc(data["branchdc"]["2"],hoa["candidate_lines"]["dc"]["23"],data["gen"]["1"]["mbase"])

    data["convdc"]["1"]=convdc(data["convdc"]["1"],hoa["nodes"]["1"],data["gen"]["1"]["mbase"])
    data["convdc"]["2"]=convdc(data["convdc"]["2"],hoa["nodes"]["2"],data["gen"]["1"]["mbase"])
    data["convdc_ne"]["1"]=convdc_ne(data["convdc_ne"]["1"],hoa["nodes"]["3"],data["gen"]["1"]["mbase"])

    #data["gen"]["1"]=gen(data["gen"]["1"],hoa["nodes"]["1"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["2"]=gen(data["gen"]["2"],-1*hoa["nodes"]["2"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["3"]=gen(data["gen"]["3"],hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["4"]=gen(data["gen"]["4"],-1*hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["5"]=gen(data["gen"]["5"],-1*hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    return data
end


function pm_dict(data,hoa)
    data["ne_branch"]["1"]=ne_branch(data["ne_branch"]["1"], hoa["candidate_lines"]["ac"]["13"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_xfo(data["ne_branch"]["1"],hoa["nodes"]["1"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_xfo(data["ne_branch"]["1"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["1"]=branch_plat(data["ne_branch"]["1"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["ne_branch"]["2"]=ne_branch(data["ne_branch"]["2"], hoa["candidate_lines"]["ac"]["23"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_xfo(data["ne_branch"]["2"],hoa["nodes"]["2"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_xfo(data["ne_branch"]["2"],hoa["nodes"]["3"]["ac"]["xfm"], data["gen"]["1"]["mbase"])
    data["ne_branch"]["2"]=branch_plat(data["ne_branch"]["2"], hoa["nodes"]["3"]["ac"]["plat"], data["gen"]["1"]["mbase"])

    data["branchdc_ne"]["1"]=branchdc_ne(data["branchdc_ne"]["1"],hoa["candidate_lines"]["dc"]["12"],data["gen"]["1"]["mbase"])
    data["branchdc_ne"]["2"]=branchdc_ne(data["branchdc_ne"]["2"],hoa["candidate_lines"]["dc"]["13"],data["gen"]["1"]["mbase"])
    data["branchdc_ne"]["3"]=branchdc_ne(data["branchdc_ne"]["3"],hoa["candidate_lines"]["dc"]["23"],data["gen"]["1"]["mbase"])

    data["convdc_ne"]["1"]=convdc_ne(data["convdc_ne"]["1"],hoa["nodes"]["1"],data["gen"]["1"]["mbase"])
    data["convdc_ne"]["2"]=convdc_ne(data["convdc_ne"]["2"],hoa["nodes"]["2"],data["gen"]["1"]["mbase"])
    data["convdc_ne"]["3"]=convdc_ne(data["convdc_ne"]["3"],hoa["nodes"]["3"],data["gen"]["1"]["mbase"])

    data["gen"]["1"]=gen(data["gen"]["1"],hoa["nodes"]["1"]["mva"]/data["gen"]["1"]["mbase"])
    data["gen"]["2"]=gen(data["gen"]["2"],-1*hoa["nodes"]["2"]["mva"]/data["gen"]["1"]["mbase"])
    data["gen"]["3"]=gen(data["gen"]["3"],hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["4"]=gen(data["gen"]["4"],-1*hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    #data["gen"]["5"]=gen(data["gen"]["5"],-1*hoa["nodes"]["3"]["mva"]/data["gen"]["1"]["mbase"])
    return data
end

function gen(data,mva)
    data["pmax"]=mva
    data["pmin"]=mva
    return data
end

function pm_dict_IConly(data,hoa)
    data["branchdc_ne"]["1"]=branchdc_ne(data["branchdc_ne"]["1"],hoa["candidate_lines"]["dc"]["12"],data["gen"]["1"]["mbase"])
    data["convdc_ne"]["1"]=convdc_ne(data["convdc_ne"]["1"],hoa["nodes"]["1"],data["gen"]["1"]["mbase"])
    data["convdc_ne"]["2"]=convdc_ne(data["convdc_ne"]["2"],hoa["nodes"]["2"],data["gen"]["1"]["mbase"])
    return data
end

function converterdc_ne(conv,nde,mbase)
    #NOTE change
    conv["cost"]=nde["xfm_off"].costs.cpx_p+nde["xfm_off"].costs.cpx_i+nde["conv_off"].costs.cpx+nde["plat"].costs.cpx+nde["cableWF"].costs.cpx_p+nde["cableWF"].costs.cpx_i;
    conv["Pacmax"]=nde["conv_off"].mva/mbase
    conv["Pacmin"]=-1*nde["conv_off"].mva/mbase
    conv["Qacmax"]=nde["conv_off"].mva/mbase
    conv["Qacmin"]=-1*nde["conv_off"].mva/mbase
    return conv
end

function convdc_ne(conv,nde,mbase)
    #NOTE change
    conv["cost"]=nde["dc"]["xfm"].costs.cpx_p+nde["dc"]["xfm"].costs.cpx_i+nde["dc"]["conv"].costs.cpx
    #conv["cost"]=0
    if (haskey(nde["dc"],"plat")); conv["cost"]=conv["cost"]+nde["dc"]["plat"].costs.cpx;end
    #if (haskey(nde["dc"],"plat")); conv["cost"]=0;end
    #########
    conv["Pacmax"]=nde["mva"]/mbase
    conv["Pacmin"]=-1*nde["mva"]/mbase
    conv["Qacmax"]=nde["mva"]/mbase
    conv["Qacmin"]=-1*nde["mva"]/mbase
    return conv
end

function branchdc_off(brc,dcoss,mbase)
    brc["cost"]=brc["cost"]+dcoss["xfm_off"].costs.cpx_p+dcoss["xfm_off"].costs.cpx_i+dcoss["conv_off"].costs.cpx+dcoss["plat"].costs.cpx
    return brc
end

function branchdc_on(brc,dcss,mbase)
    brc["cost"]=brc["cost"]+dcss["xfm_on"].costs.cpx_p+dcss["xfm_on"].costs.cpx_i+dcss["conv_on"].costs.cpx
    return brc
end

function convdc(conv,nde,mbase)
    conv["Pacmax"]=nde["mva"]/mbase
    conv["Pacmin"]=-1*nde["mva"]/mbase
    conv["Qacmax"]=nde["mva"]/mbase
    conv["Qacmin"]=-1*nde["mva"]/mbase
    return conv
end

function converterdc(conv,nde,mbase)
    conv["Pacmax"]=nde.mva/mbase
    conv["Pacmin"]=-1*nde.mva/mbase
    conv["Qacmax"]=nde.mva/mbase
    conv["Qacmin"]=-1*nde.mva/mbase
    return conv
end

function branchdc_ne(b_ne,dc_can,mbase)
    b_ne["rateA"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateB"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateC"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["r"]=dc_can.elec.ohm

    #NOTE change
    b_ne["cost"]=dc_can.costs.cpx_p+dc_can.costs.cpx_i
    #b_ne["cost"]=0
    #######
    return b_ne
end

function branchdc(b_ne,dc_can,mbase)
    b_ne["rateA"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateB"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateC"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["r"]=dc_can.elec.ohm
    return b_ne
end

function ne_branch(ne_b,ac_can,mbase)
    ne_b["mva"]=ac_can.num*ac_can.elec.mva
    ne_b["rate_a"]=ac_can.num*ac_can.elec.mva/mbase
    ne_b["rate_b"]=ac_can.num*ac_can.elec.mva/mbase
    ne_b["rate_c"]=ac_can.num*ac_can.elec.mva/mbase
    ne_b["br_r"]=ac_can.elec.ohm
    ne_b["br_x"]=ac_can.elec.xl
    #data["ne_branch"][datai]["construction_cost"]=hoa["candidate_lines"]["ac"][hoai].costs.cpx_p+hoa["candidate_lines"]["ac"][hoai].costs.cpx_i+hoa["candidate_lines"]["ac"][hoai].costs.qc+hoa["candidate_lines"]["ac"][hoai].costs.sg
    #ne_b["construction_cost"]=(ac_can.costs.cpx_p+ac_can.costs.cpx_i)/mbase

    #NOTE change
    ne_b["construction_cost"]=ac_can.costs.cpx_p+ac_can.costs.cpx_i
    #ne_b["construction_cost"]=0
    ##########
    return ne_b
end

function branch_xfo(ne_b,ne_xfo,mbase)
    #ne_b["construction_cost"]=ne_b["construction_cost"]+(ne_xfo.costs.cpx_p+ne_xfo.costs.cpx_i)/mbase

    #NOTE change
    ne_b["construction_cost"]=ne_b["construction_cost"]+ne_xfo.costs.cpx_p+ne_xfo.costs.cpx_i
    #ne_b["construction_cost"]=ne_b["construction_cost"]+0
    ##########
    return ne_b
end

function branch_plat(ne_b,ne_plat,mbase)
    #ne_b["construction_cost"]=ne_b["construction_cost"]+(ne_plat.costs.cpx)/mbase

    #NOTE change
    ne_b["construction_cost"]=ne_b["construction_cost"]+ne_plat.costs.cpx
    #ne_b["construction_cost"]=ne_b["construction_cost"]+0
    ##########
    return ne_b
end
#calculates candidate equipment and sorts in struct
function hoa_datastruct(ic_mva,owpp_mva,ic_length,owpp_km)
    #location
    uk_pcc=xy();uk_pcc.x=ic_length;uk_pcc.y=0;
    be_pcc=xy();be_pcc.x=0;be_pcc.y=0;
    owpp_xy=xy();owpp_xy.x=owpp_km;owpp_xy.y=0
    #datastructure
    hoa=Dict{String,Any}();
    push!(hoa,"nodes"=>Dict{String,Any}());

    [push!(hoa["nodes"],n=>Dict{String,Any}()) for n in ["1","2","3"]]
    #position of nodes
    push!(hoa["nodes"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["2"],"xy"=>uk_pcc);
    push!(hoa["nodes"]["3"],"xy"=>owpp_xy);
    #rating of nodes
    push!(hoa["nodes"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["2"],"mva"=>ic_mva);
    push!(hoa["nodes"]["3"],"mva"=>owpp_mva);

    [push!(hoa,"ac"=>Dict{String,Any}()) for hoa in [hoa["nodes"]["1"],hoa["nodes"]["2"],hoa["nodes"]["3"]]]
    [push!(hoa,"dc"=>Dict{String,Any}()) for hoa in [hoa["nodes"]["1"],hoa["nodes"]["2"],hoa["nodes"]["3"]]]
    #node 1 - ac
    push!(hoa["nodes"]["1"]["ac"],"xfm"=>on_xfm(hoa["nodes"]["3"]["mva"]))
    #node 1 - dc
    push!(hoa["nodes"]["1"]["dc"],"xfm"=>on_xfm(hoa["nodes"]["1"]["mva"]))
    push!(hoa["nodes"]["1"]["dc"],"conv"=>on_conv(hoa["nodes"]["1"]["mva"]))
    #node 2 - ac
    push!(hoa["nodes"]["2"]["ac"],"xfm"=>on_xfm(hoa["nodes"]["3"]["mva"]))
    #node 2 - dc
    push!(hoa["nodes"]["2"]["dc"],"xfm"=>on_xfm(hoa["nodes"]["2"]["mva"]))
    push!(hoa["nodes"]["2"]["dc"],"conv"=>on_conv(hoa["nodes"]["2"]["mva"]))
    #node 3 -ac
    push!(hoa["nodes"]["3"]["ac"],"plat"=>AC_plat(hoa["nodes"]["3"]["mva"]))
    push!(hoa["nodes"]["3"]["ac"],"xfm"=>off_xfm(hoa["nodes"]["3"]["mva"]))
    #node 3 -dc
    push!(hoa["nodes"]["3"]["dc"],"plat"=>DC_plat(hoa["nodes"]["3"]["mva"]))
    push!(hoa["nodes"]["3"]["dc"],"xfm"=>off_xfm(hoa["nodes"]["3"]["mva"]))
    push!(hoa["nodes"]["3"]["dc"],"conv"=>off_conv(hoa["nodes"]["3"]["mva"]))



    #lines
    push!(hoa,"candidate_lines"=>Dict{String,Any}());
    push!(hoa["candidate_lines"],"ac"=>Dict{String,Any}());
    push!(hoa["candidate_lines"],"dc"=>Dict{String,Any}());
    #candidate ac lines
    push!(hoa["candidate_lines"]["ac"],"13"=>AC_cbl(owpp_mva, owpp_km));
    push!(hoa["candidate_lines"]["ac"],"23"=>AC_cbl(owpp_mva, ic_length-owpp_km));
    #candidate dc lines
    push!(hoa["candidate_lines"]["dc"],"13"=>DC_cbl(ic_mva, owpp_km));
    push!(hoa["candidate_lines"]["dc"],"23"=>DC_cbl(ic_mva, ic_length-owpp_km));
    push!(hoa["candidate_lines"]["dc"],"12"=>DC_cbl(ic_mva, ic_length));
    return hoa
end

function pm_dict_candidateICs(data,hoa,candidates)
    #WF AC candidate lines
    for (i,k) in enumerate(vcat(candidates[2],candidates[2]))
        if (i<=length(candidates[2]))
            data["ne_branch"][string(i)]=ne_branch(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["cableBE"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_on"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_off"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_plat(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["plat"], data["gen"]["1"]["mbase"])
        else
            data["ne_branch"][string(i)]=ne_branch(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["cableUK"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_on"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_off"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_plat(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["plat"], data["gen"]["1"]["mbase"])
        end
    end
    #WF DC
    #=for (i,k) in enumerate(candidates[2])
        data["branchdc_ne"][string(i)]=branchdc_ne(data["branchdc_ne"][string(i)],hoa["candidates"]["wf"]["dc"][string(k)]["cableWF"],data["gen"]["1"]["mbase"])
        data["branchdc_ne"][string(i)]=branchdc_off(data["branchdc_ne"][string(i)],hoa["candidates"]["wf"]["dc"][string(k)],data["gen"]["1"]["mbase"])
    end=#

    for (i,k) in enumerate(vcat(candidates[1],candidates[1]))
        if (i<=length(candidates[2]))
            data["branchdc_ne"][string(i)]=branchdc_ne(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)]["cableBE"],data["gen"]["1"]["mbase"])
            data["branchdc_ne"][string(i)]=branchdc_on(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)],data["gen"]["1"]["mbase"])
        else
            data["branchdc_ne"][string(i)]=branchdc_ne(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)]["cableUK"],data["gen"]["1"]["mbase"])
            data["branchdc_ne"][string(i)]=branchdc_on(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)],data["gen"]["1"]["mbase"])
        end
    end

    for (i,k) in enumerate(candidates[2])
        data["convdc_ne"][string(i)]=converterdc_ne(data["convdc_ne"][string(i)],hoa["candidates"]["wf"]["dc"][string(k)],data["gen"]["1"]["mbase"])
    end
    size_set=[];for (k,v) in hoa["candidates"]["ic"]; push!(size_set,parse(Float64,k));end
    mx_key=string(maximum(size_set))
    data["convdc"]["1"]=converterdc(data["convdc"]["1"],hoa["candidates"]["ic"][mx_key]["conv_on"],data["gen"]["1"]["mbase"])
    data["convdc"]["2"]=converterdc(data["convdc"]["2"],hoa["candidates"]["ic"][mx_key]["conv_on"],data["gen"]["1"]["mbase"])


    #=data["gen"]["1"]=gen(data["gen"]["1"],ic_mva/data["gen"]["1"]["mbase"])
    data["gen"]["2"]=gen(data["gen"]["2"],-1*ic_mva/data["gen"]["1"]["mbase"])
    data["gen"]["3"]=gen(data["gen"]["3"],ic_mva/data["gen"]["1"]["mbase"])
    data["gen"]["4"]=gen(data["gen"]["4"],-1*ic_mva/data["gen"]["1"]["mbase"])
    data["gen"]["5"]=gen(data["gen"]["5"],owpp_mva/data["gen"]["1"]["mbase"])=#

    return data
end

function additional_candidates(data,candidates)
    #AC
    c2be=deepcopy(data["ne_branch"]["1"]);c2uk=deepcopy(data["ne_branch"]["2"])
    data["ne_branch"]=Dict{String,Any}()
    #WF
    for i=1:1:length(candidates[2]); c2be["source_id"][2]=i;c2be["index"]=i; push!(data["ne_branch"],string(i)=>deepcopy(c2be)); end
    for i=length(candidates[2])+1:1:length(candidates[2])*2; c2uk["source_id"][2]=i;c2uk["index"]=i; push!(data["ne_branch"],string(i)=>deepcopy(c2uk)); end
    #DC
    #WF
    convwf=deepcopy(data["convdc_ne"]["1"])
    data["convdc_ne"]=Dict{String,Any}()
    for i=1:1:length(candidates[2]); convwf["source_id"][2]=i;convwf["index"]=i; push!(data["convdc_ne"],string(i)=>deepcopy(convwf)); end
    #IC
    c2be=deepcopy(data["branchdc_ne"]["1"]);c2uk=deepcopy(data["branchdc_ne"]["2"])
    data["branchdc_ne"]=Dict{String,Any}()
    for i=1:1:length(candidates[1]); c2be["source_id"][2]=i;c2be["index"]=i; push!(data["branchdc_ne"],string(i)=>deepcopy(c2be)); end
    for i=length(candidates[1])+1:1:length(candidates[1])*2; c2uk["source_id"][2]=i;c2uk["index"]=i; push!(data["branchdc_ne"],string(i)=>deepcopy(c2uk)); end
    return data
end

#calculates candidate equipment and sorts in struct
function hoa_datastruct_candidateICs(ic_mva,owpp_mva,ic_length,owpp_km,candidates)
    #location
    uk_pcc=xy();uk_pcc.x=ic_length;uk_pcc.y=0;
    be_pcc=xy();be_pcc.x=0;be_pcc.y=0;
    owpp_xy=xy();owpp_xy.x=owpp_km;owpp_xy.y=0
    #datastructure
    hoa=Dict{String,Any}();
    push!(hoa,"nodes"=>Dict{String,Any}());
    push!(hoa["nodes"],"ac"=>Dict{String,Any}());
    push!(hoa["nodes"],"dc"=>Dict{String,Any}());

    [push!(hoa["nodes"]["ac"],n=>Dict{String,Any}()) for n in ["1","2","3"]]
    [push!(hoa["nodes"]["dc"],n=>Dict{String,Any}()) for n in ["1","2","3","4"]]
    #position of nodes
    push!(hoa["nodes"]["ac"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["ac"]["2"],"xy"=>uk_pcc);
    push!(hoa["nodes"]["ac"]["3"],"xy"=>owpp_xy);
    #position of nodes
    push!(hoa["nodes"]["dc"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["dc"]["2"],"xy"=>uk_pcc);
    push!(hoa["nodes"]["dc"]["3"],"xy"=>owpp_xy);
    push!(hoa["nodes"]["dc"]["4"],"xy"=>owpp_xy);
    #rating of nodes
    push!(hoa["nodes"]["ac"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["ac"]["2"],"mva"=>ic_mva);
    push!(hoa["nodes"]["ac"]["3"],"mva"=>owpp_mva);
    #rating of nodes
    push!(hoa["nodes"]["dc"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["dc"]["2"],"mva"=>ic_mva);
    push!(hoa["nodes"]["dc"]["3"],"mva"=>owpp_mva);
    push!(hoa["nodes"]["dc"]["4"],"mva"=>owpp_mva);
    #lines
    push!(hoa,"candidates"=>Dict{String,Any}());
    push!(hoa["candidates"],"wf"=>Dict{String,Any}());
    push!(hoa["candidates"],"ic"=>Dict{String,Any}());
    push!(hoa["candidates"]["wf"],"ac"=>Dict{String,Any}());
    push!(hoa["candidates"]["wf"],"dc"=>Dict{String,Any}());

    [push!(hoa["candidates"]["wf"]["ac"],string(p)=>Dict{String,Any}()) for p in candidates[2]]
    [push!(hoa["candidates"]["wf"]["dc"],string(p)=>Dict{String,Any}()) for p in candidates[2]]
    [push!(hoa["candidates"]["ic"],string(p)=>Dict{String,Any}()) for p in candidates[1]]

    #wf-ac
    [push!(hoa["candidates"]["wf"]["ac"][k],"xfm_on"=>on_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"xfm_off"=>off_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"plat"=>AC_plat(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"xfm"=>on_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"cableBE"=>AC_cbl(owpp_mva*parse(Float64,k), owpp_km)) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"cableUK"=>AC_cbl(owpp_mva*parse(Float64,k), ic_length-owpp_km)) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    #wf-dc
    [push!(hoa["candidates"]["wf"]["dc"][k],"xfm_off"=>off_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    [push!(hoa["candidates"]["wf"]["dc"][k],"conv_off"=>off_conv(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    [push!(hoa["candidates"]["wf"]["dc"][k],"plat"=>DC_plat(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    [push!(hoa["candidates"]["wf"]["dc"][k],"cableWF"=>DC_cbl(owpp_mva/2*parse(Float64,k), 1)) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    #ic dc
    [push!(hoa["candidates"]["ic"][k],"xfm_on"=>on_xfm(ic_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["ic"]]
    [push!(hoa["candidates"]["ic"][k],"conv_on"=>on_conv(ic_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["ic"]]
    [push!(hoa["candidates"]["ic"][k],"cableBE"=>DC_cbl(ic_mva/2*parse(Float64,k), owpp_km)) for (k,p) in hoa["candidates"]["ic"]]
    [push!(hoa["candidates"]["ic"][k],"cableUK"=>DC_cbl(ic_mva/2*parse(Float64,k), ic_length-owpp_km)) for (k,p) in hoa["candidates"]["ic"]]

    return hoa
end

function hoa_datastruct_owpp_only(ic_mva,owpp_mva,ic_length,owpp_km)
    #location
    uk_pcc=xy();uk_pcc.x=ic_length;uk_pcc.y=0;
    be_pcc=xy();be_pcc.x=0;be_pcc.y=0;
    owpp_xy=xy();owpp_xy.x=owpp_km;owpp_xy.y=0
    #datastructure
    hoa=Dict{String,Any}();
    push!(hoa,"nodes"=>Dict{String,Any}());

    [push!(hoa["nodes"],n=>Dict{String,Any}()) for n in ["1","2","3"]]
    #position of nodes
    push!(hoa["nodes"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["2"],"xy"=>uk_pcc);
    push!(hoa["nodes"]["3"],"xy"=>owpp_xy);
    #rating of nodes
    push!(hoa["nodes"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["2"],"mva"=>ic_mva);
    push!(hoa["nodes"]["3"],"mva"=>owpp_mva);

    [push!(hoa,"ac"=>Dict{String,Any}()) for hoa in [hoa["nodes"]["1"],hoa["nodes"]["2"],hoa["nodes"]["3"]]]
    [push!(hoa,"dc"=>Dict{String,Any}()) for hoa in [hoa["nodes"]["1"],hoa["nodes"]["2"],hoa["nodes"]["3"]]]
    #node 1 - ac
    push!(hoa["nodes"]["1"]["ac"],"xfm"=>on_xfm(hoa["nodes"]["3"]["mva"]))
    #node 1 - dc
    push!(hoa["nodes"]["1"]["dc"],"xfm"=>on_xfm(hoa["nodes"]["1"]["mva"]))
    push!(hoa["nodes"]["1"]["dc"],"conv"=>on_conv(hoa["nodes"]["1"]["mva"]))
    #node 2 - ac
    push!(hoa["nodes"]["2"]["ac"],"xfm"=>on_xfm(hoa["nodes"]["3"]["mva"]))
    #node 2 - dc
    push!(hoa["nodes"]["2"]["dc"],"xfm"=>on_xfm(hoa["nodes"]["2"]["mva"]))
    push!(hoa["nodes"]["2"]["dc"],"conv"=>on_conv(hoa["nodes"]["2"]["mva"]))
    #node 3 -ac
    push!(hoa["nodes"]["3"]["ac"],"plat"=>AC_plat(hoa["nodes"]["3"]["mva"]))
    push!(hoa["nodes"]["3"]["ac"],"xfm"=>off_xfm(hoa["nodes"]["3"]["mva"]))
    #node 3 -dc
    push!(hoa["nodes"]["3"]["dc"],"plat"=>DC_plat(hoa["nodes"]["3"]["mva"]))
    push!(hoa["nodes"]["3"]["dc"],"xfm"=>off_xfm(hoa["nodes"]["3"]["mva"]))
    push!(hoa["nodes"]["3"]["dc"],"conv"=>off_conv(hoa["nodes"]["3"]["mva"]))



    #lines
    push!(hoa,"candidate_lines"=>Dict{String,Any}());
    push!(hoa["candidate_lines"],"ac"=>Dict{String,Any}());
    push!(hoa["candidate_lines"],"dc"=>Dict{String,Any}());
    #candidate ac lines
    push!(hoa["candidate_lines"]["ac"],"13"=>[AC_cbl(owpp_mva, owpp_km),AC_cbl(owpp_mva*0.95, owpp_km),AC_cbl(owpp_mva*0.90, owpp_km),AC_cbl(owpp_mva*0.85, owpp_km),AC_cbl(owpp_mva*0.80, owpp_km)]);
    push!(hoa["candidate_lines"]["ac"],"23"=>AC_cbl(owpp_mva, ic_length-owpp_km));
    #candidate dc lines
    push!(hoa["candidate_lines"]["dc"],"13"=>DC_cbl(ic_mva, owpp_km));
    push!(hoa["candidate_lines"]["dc"],"23"=>DC_cbl(ic_mva, ic_length-owpp_km));
    push!(hoa["candidate_lines"]["dc"],"12"=>DC_cbl(ic_mva, ic_length));
    return hoa
end

function hoa_datastruct_IConly(ic_mva,owpp_mva,ic_length,owpp_km)
    #location
    uk_pcc=xy();uk_pcc.x=ic_length;uk_pcc.y=0;
    be_pcc=xy();be_pcc.x=0;be_pcc.y=0;

    #datastructure
    hoa=Dict{String,Any}();
    push!(hoa,"nodes"=>Dict{String,Any}());

    [push!(hoa["nodes"],n=>Dict{String,Any}()) for n in ["1","2"]]
    #position of nodes
    push!(hoa["nodes"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["2"],"xy"=>uk_pcc);

    #rating of nodes
    push!(hoa["nodes"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["2"],"mva"=>ic_mva);

    [push!(hoa,"dc"=>Dict{String,Any}()) for hoa in [hoa["nodes"]["1"],hoa["nodes"]["2"]]]

    #node 1 - dc
    push!(hoa["nodes"]["1"]["dc"],"xfm"=>on_xfm(hoa["nodes"]["1"]["mva"]))
    push!(hoa["nodes"]["1"]["dc"],"conv"=>on_conv(hoa["nodes"]["1"]["mva"]))

    #node 2 - dc
    push!(hoa["nodes"]["2"]["dc"],"xfm"=>on_xfm(hoa["nodes"]["2"]["mva"]))
    push!(hoa["nodes"]["2"]["dc"],"conv"=>on_conv(hoa["nodes"]["2"]["mva"]))



    #lines
    push!(hoa,"candidate_lines"=>Dict{String,Any}());
    push!(hoa["candidate_lines"],"dc"=>Dict{String,Any}());

    #candidate dc lines
    push!(hoa["candidate_lines"]["dc"],"12"=>DC_cbl(ic_mva, ic_length));
    return hoa
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

function rand_set_days(ukbe_data,dn,wn,mn)
    sampled_set=DataFrame()
    selected_months=[]; for j=1:mn; months=[i for i=1:12 if !issubset(i,selected_months)]; push!(selected_months,rand(months));end;sort!(selected_months)
    for mnth in selected_months;#println("month: "*string(mnth))
        hours_in_month=filter(e->month(e.time_stamp)==mnth,ukbe_data)
        days_in_month=round(Int64,length(hours_in_month[!,:time_stamp])/24)
        weeks_in_month=ceil(Int64,days_in_month/7);selected_weeks=[];
        #for w=1:wn; weeks=[i for i=1:weeks_in_month if !issubset(i,selected_weeks)]; push!(selected_weeks,rand(weeks));end;sort!(selected_weeks);selected_days=[];
        selected_weeks=wn;selected_days=[];
        for wek in selected_weeks;#println(wek)
            for d=1:dn; days=[i for i in (wek-1)*7+1:min((wek)*7,days_in_month) if !issubset(i,selected_days)]; if (length(days)>0); push!(selected_days,rand(days));end;end;sort!(selected_days);end
        for dy in selected_days;
            hours_in_day=filter(e->day(e.time_stamp)==dy,hours_in_month)
            sampled_set=vcat(sampled_set,hours_in_day)
        end
    end
    sort!(sampled_set)
    return sampled_set
end


function get_profile_data_sets(d1,d2,data, n, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];

    df0=DataFrame("time_stamp"=>[],"daprice"=>[],"idprice"=>[],"wind"=>[],"regup_price"=>[],"regup_mwh"=>[],"regdwn_price"=>[],"regdwn_mwh"=>[])
    df1=DataFrame("time_stamp"=>[],"daprice"=>[],"idprice"=>[],"wind"=>[],"regup_price"=>[],"regup_mwh"=>[],"regdwn_price"=>[],"regdwn_mwh"=>[])#MWh_up,EUR_up,MWh_dwn,EUR_dwn
    z0_data=CSV.read("./test/data/cordoba/input/"*d1*".csv", DataFrames.DataFrame)
    z1_data=CSV.read("./test/data/cordoba/input/"*d2*".csv", DataFrames.DataFrame)
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


function get_profile_data_UKBE(data,dn,wn,mn, scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    if haskey(scenario, "mc")
        monte_carlo = scenario["mc"]
    else
        monte_carlo = false
    end

    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]
        ukbe_data=CSV.read("./test/data/cordoba/input/ukbe_ts_2.csv", DataFrames.DataFrame)
        #ukbe_data=ukbe_data[shuffle(1:nrow(ukbe_data)),scenario["hours"]),:]
        ############# random n hours ####################
        #ukbe_data=ukbe_data[shuffle(1:nrow(ukbe_data))[1:scenario["hours"]], :]
        #sort!(ukbe_data)

        ############### random maintaining days specify dn,wn,mn
        ukbe_data=rand_set_days(ukbe_data,dn,wn,mn)
        #=it=0
        tol=0.001
        while (it<1000 && !(sum(ukbe_data.be_wind)/length(ukbe_data.be_wind)+tol>0.362 && sum(ukbe_data.be_wind)/length(ukbe_data.be_wind)-tol<0.362 && std(ukbe_data.be_wind)+tol>0.308 && std(ukbe_data.be_wind)-tol<0.308))
            ukbe_data=CSV.read("./test/data/cordoba/input/ukbe_ts.csv", DataFrames.DataFrame)
            #ukbe_data=rand_set_days(ukbe_data,dn,wn,mn)
            ukbe_data=ukbe_data[shuffle(1:nrow(ukbe_data))[1:scenario["hours"]], :]
            sort!(ukbe_data)
            it=it+1
        end
        if it<1000
            println("cf: "*string(sum(ukbe_data.be_wind)/length(ukbe_data.be_wind)))
            println("std: "*string(std(ukbe_data.be_wind)))
        else
            println("not found")
        end=#
        #ukbe_data.time_stamp = DateTime.(ukbe_data.time_stamp,dateformat"dd.mm.yyyy HH:MM")
        #ukbe_data.time_stamp = DateTime.(ukbe_data.time_stamp,dateformat"yyyy-mm-ddTHH:MM:SS")
        if monte_carlo == false
            utch_idx = Dates.epochms2datetime(scnr["start"])
            #ukbe_range=ukbe_data[(ukbe_data[:time_stamp].>=utch_idx).&(ukbe_data[:time_stamp].<utch_idx+Hour(scenario["hours"])),:]
            #ukbe_range=ukbe_data[(ukbe_data[!,:time_stamp].>=utch_idx).&(ukbe_data[!,:time_stamp].<utch_idx+Hour(scenario["hours"])),:]
            ukbe_range=ukbe_data
            #time_stamp,be_costid,be_eumwh,be_wind,uk_costid,uk_eumwh,uk_wind,losses
            #windgenprofile_beuk=[ukbe_range[:be_wind]';ukbe_range[:uk_wind]']
            #gencostid_beuk=[ukbe_range[:be_costid]';ukbe_range[:uk_costid]']
            #gencost_beuk=[ukbe_range[:be_eumwh]';ukbe_range[:uk_eumwh]']
            #ic_nflow_losses=[ukbe_range[:net_flows]';ukbe_range[:losses]']
            windgenprofile_beuk=[ukbe_range[!,:be_wind]';ukbe_range[!,:uk_wind]']
            gencostid_beuk=[ukbe_range[!,:be_costid]';ukbe_range[!,:uk_costid]']
            gencost_beuk=[ukbe_range[!,:be_eumwh]';ukbe_range[!,:uk_eumwh]']
            ic_nflow_losses=[ukbe_range[!,:net_flows]';ukbe_range[!,:losses]']
        else

        end
        scenario["hours"]=length(ic_nflow_losses[1,:])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, ic_nflow_losses
end

function get_profile_data_UKBE(data,scenario = Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
    if haskey(scenario, "mc")
        monte_carlo = scenario["mc"]
    else
        monte_carlo = false
    end

    for (s, scnr) in scenario["sc_years"]
        start_idx = (parse(Int, s) - 1) * scenario["hours"]
        ukbe_data=CSV.read("./test/data/cordoba/input/ukbe_ts_2.csv", DataFrames.DataFrame)
        #ukbe_data=ukbe_data[shuffle(1:nrow(ukbe_data)),scenario["hours"]),:]
        ############# random n hours ####################
        #ukbe_data=ukbe_data[shuffle(1:nrow(ukbe_data))[1:scenario["hours"]], :]
        sort!(ukbe_data)

        ############### random maintaining days specify dn,wn,mn
        #ukbe_data=rand_set_days(ukbe_data,dn,wn,mn)
        #ukbe_data.time_stamp = DateTime.(ukbe_data.time_stamp,dateformat"dd.mm.yyyy HH:MM")
        #ukbe_data.time_stamp = DateTime.(ukbe_data.time_stamp,dateformat"yyyy-mm-ddTHH:MM:SS")
        if monte_carlo == false
            utch_idx = Dates.epochms2datetime(scnr["start"])
            #ukbe_range=ukbe_data[(ukbe_data[:time_stamp].>=utch_idx).&(ukbe_data[:time_stamp].<utch_idx+Hour(scenario["hours"])),:]
            ukbe_range=ukbe_data[(ukbe_data[!,:time_stamp].>=utch_idx).&(ukbe_data[!,:time_stamp].<utch_idx+Hour(scenario["hours"])),:]
            #ukbe_range=ukbe_data
            #time_stamp,be_costid,be_eumwh,be_wind,uk_costid,uk_eumwh,uk_wind,losses
            #windgenprofile_beuk=[ukbe_range[:be_wind]';ukbe_range[:uk_wind]']
            #gencostid_beuk=[ukbe_range[:be_costid]';ukbe_range[:uk_costid]']
            #gencost_beuk=[ukbe_range[:be_eumwh]';ukbe_range[:uk_eumwh]']
            #ic_nflow_losses=[ukbe_range[:net_flows]';ukbe_range[:losses]']
            windgenprofile_beuk=[ukbe_range[!,:be_wind]';ukbe_range[!,:uk_wind]']
            gencostid_beuk=[ukbe_range[!,:be_costid]';ukbe_range[!,:uk_costid]']
            gencost_beuk=[ukbe_range[!,:be_eumwh]';ukbe_range[!,:uk_eumwh]']
            ic_nflow_losses=[ukbe_range[!,:net_flows]';ukbe_range[!,:losses]']
        else

        end
        scenario["hours"]=length(ic_nflow_losses[1,:])
        data["scenario"][s] = Dict()
        data["scenario_prob"][s] = scnr["probability"]
        for h in 1 : scenario["hours"]
            network = start_idx + h
            data["scenario"][s]["$h"] = network
        end

    end
    # Return info
    return data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, ic_nflow_losses
end



function create_profile_UKBE(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses)
    #e2me=1000000/data["baseMVA"]#into ME/PU
    e2me=1#into ME/PU
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
        if (nFlow_losses[1, d]>=0)#net flow is from BE to UK
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]

            #extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["1"]["pmin"][1, d] = 0
            ##############

            #extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]]
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]
        #load generator

            extradata["gen"]["2"]["pmax"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            #extradata["gen"]["2"]["pmax"][1, d] = 0
            ##############

            #extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]#NOTE HERE
            extradata["gen"]["2"]["pmin"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            ##############

            #extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d] * -1]
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
        else#net flow is from UK to BE
            #load generator

            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            #extradata["gen"]["1"]["pmax"][1, d] = 0
            ##############

            #extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]#NOTE HERE
            extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            ##############

            #extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]
        #source generator
            extradata["gen"]["2"]["pmax"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]

            #extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["2"]["pmin"][1, d] = 0
            ##############

            #extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]]
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
        end
        extradata["gen"]["3"]["pmax"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]

        #extradata["gen"]["3"]["pmin"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        extradata["gen"]["3"]["pmin"][1, d] = 0
        ##############

        #extradata["gen"]["3"]["cost"][d] = [0.0]
        wnd_cost=min(0.0,gencost_beuk[1,d],gencost_beuk[2,d])
        #wnd_cost=min(0.0,21.0,1.0)
        extradata["gen"]["3"]["cost"][d] = [wnd_cost/e2me,0.0]

        extradata["gen"]["4"]["pmax"][1, d] = 0.0
        #NOTE changed
        ##############
        extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        #extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d] * abs(nFlow_losses[1, d])
        #extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
        extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]

        extradata["gen"]["5"]["pmax"][1, d] = 0.0
        #NOTE changed
        ##############
        extradata["gen"]["5"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        #extradata["gen"]["5"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d] * abs(nFlow_losses[1, d])
        #extradata["gen"]["5"]["cost"][d] = [gencost_beuk[2,d] * -1.0]
        extradata["gen"]["5"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
    end
    return extradata
end

function create_profile_UKBE_0835(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses, ic_mva, owpp_mva)
    e2me=1000000/data["baseMVA"]#into ME/PU
    da=1;id=0
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
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/data["baseMVA"]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/data["baseMVA"]
            extradata["gen"]["2"]["cost"][d] = [da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,0]

        #Intra ahead BE
        #source generator
            #=extradata["gen"]["3"]["pmax"][1, d] = id*ic_mva/data["baseMVA"]
            extradata["gen"]["3"]["pmin"][1, d] = 0
            extradata["gen"]["3"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]
        #load generator
            extradata["gen"]["4"]["pmax"][1, d] = -id*ic_mva/data["baseMVA"]
            extradata["gen"]["4"]["pmin"][1, d] = -id*ic_mva/data["baseMVA"]
            extradata["gen"]["4"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]=#


            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/data["baseMVA"]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me,0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/data["baseMVA"]
                extradata["gen"]["4"]["cost"][d] = [da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me,0]

            #Intra ahead UK
            #source generator
                #=extradata["gen"]["7"]["pmax"][1, d] = id*ic_mva/data["baseMVA"]
                extradata["gen"]["7"]["pmin"][1, d] = 0
                extradata["gen"]["7"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]
            #load generator
                extradata["gen"]["8"]["pmax"][1, d] = -id*ic_mva/data["baseMVA"]
                extradata["gen"]["8"]["pmin"][1, d] = -id*ic_mva/data["baseMVA"]
                extradata["gen"]["8"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]=#

                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["5"]["pmin"][1, d] = 0

                #extradata["gen"]["3"]["cost"][d] = [0.0]
                #wnd_cost=min(0.0,da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me)
                #wnd_cost=min(0.0,21.0,1.0)
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]

            #Day ahead BE
            #source generator
                extradata["gen"]["6"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["6"]["pmin"][1, d] = 0
                extradata["gen"]["6"]["cost"][d] = [da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,0]
            #load generator
                extradata["gen"]["7"]["pmax"][1, d] = 0
                extradata["gen"]["7"]["pmin"][1, d] = -owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["7"]["cost"][d] = [da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,0]

            #Intra ahead BE
            #source generator
                #=extradata["gen"]["12"]["pmax"][1, d] = id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["12"]["pmin"][1, d] = 0
                extradata["gen"]["12"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]
            #load generator
                extradata["gen"]["13"]["pmax"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["13"]["pmin"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["13"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]=#


                #Day ahead UK
                #source generator
                    extradata["gen"]["8"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["8"]["pmin"][1, d] = 0
                    extradata["gen"]["8"]["cost"][d] = [da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me,0]
                #load generator
                    extradata["gen"]["9"]["pmax"][1, d] = 0
                    extradata["gen"]["9"]["pmin"][1, d] = -owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["9"]["cost"][d] = [da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me,0]

                #Intra ahead UK
                #source generator
                    #=extradata["gen"]["16"]["pmax"][1, d] = id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["16"]["pmin"][1, d] = 0
                    extradata["gen"]["16"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]
                #load generator
                    extradata["gen"]["17"]["pmax"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["17"]["pmin"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["17"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]=#
    end
    return extradata
end


function create_profile_UKBE_0835_update(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses, ic_mva, owpp_mva)
    e2me=1000000/data["baseMVA"]#into ME/PU
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
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/data["baseMVA"]+owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/data["baseMVA"]-owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
            extradata["gen"]["2"]["cost"][d] = [da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,0]

            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/data["baseMVA"]+owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me,0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/data["baseMVA"]-owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["4"]["cost"][d] = [da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me,0]
            #Wind generator
                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["5"]["pmin"][1, d] = 0
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]
    end
    return extradata
end


function create_profile_sets(number_of_hours, data, df0, df1)
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

                #=total_afrr_up=(best_up_mwh+worst_up_mwh)#/pu
                total_afrr_dwn=(best_dwn_mwh+worst_dwn_mwh)#/pu
                println(string(d)*" "*string(best_up_mwh)*" "*string(worst_up_mwh)*" "*string(best_dwn_mwh)*" "*string(worst_dwn_mwh))
                println(string(d)*" "*string(best_up_price)*" "*string(worst_up_price)*" "*string(best_dwn_price)*" "*string(worst_dwn_price))=#
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


function create_profile(number_of_hours, data, df)
    pu=data["baseMVA"]
    e2me=1000000/pu#into ME/PU
    da=0.835;id=0.165
    #e2me=1#into ME/PU
    #pu=1
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
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df.z0_wind[d]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [da*df.z0_daprice[d]/e2me+id*df.z0_idprice[d]/e2me,0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df.z0_wind[d]
            extradata["gen"]["2"]["cost"][d] = [da*df.z0_daprice[d]/e2me+id*df.z0_idprice[d]/e2me,0]

            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/pu+owpp_mva/pu*df.z0_wind[d]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [da*df.z1_daprice[d]/e2me+id*df.z1_idprice[d]/e2me,0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/pu-owpp_mva/pu*df.z0_wind[d]
                extradata["gen"]["4"]["cost"][d] = [da*df.z1_daprice[d]/e2me+id*df.z1_idprice[d]/e2me,0]
            #Wind generator
                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/pu*df.z0_wind[d]
                extradata["gen"]["5"]["pmin"][1, d] = 0
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]
    end
    return extradata
end



function create_profile_UKBE_0835_inv(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses, ic_mva, owpp_mva)
    e2me=1000000/data["baseMVA"]#into ME/PU
    da=0.085;id=0.165
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
            extradata["gen"]["1"]["pmax"][1, d] = ic_mva/data["baseMVA"]
            extradata["gen"]["1"]["pmin"][1, d] = 0
            extradata["gen"]["1"]["cost"][d] = [-1*(da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me),0]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = 0
            extradata["gen"]["2"]["pmin"][1, d] = -ic_mva/data["baseMVA"]
            extradata["gen"]["2"]["cost"][d] = [-1*(da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me),0]

        #Intra ahead BE
        #source generator
            #=extradata["gen"]["3"]["pmax"][1, d] = id*ic_mva/data["baseMVA"]
            extradata["gen"]["3"]["pmin"][1, d] = 0
            extradata["gen"]["3"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]
        #load generator
            extradata["gen"]["4"]["pmax"][1, d] = -id*ic_mva/data["baseMVA"]
            extradata["gen"]["4"]["pmin"][1, d] = -id*ic_mva/data["baseMVA"]
            extradata["gen"]["4"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]=#


            #Day ahead UK
            #source generator
                extradata["gen"]["3"]["pmax"][1, d] = ic_mva/data["baseMVA"]
                extradata["gen"]["3"]["pmin"][1, d] = 0
                extradata["gen"]["3"]["cost"][d] = [-1*(da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me),0]
            #load generator
                extradata["gen"]["4"]["pmax"][1, d] = 0
                extradata["gen"]["4"]["pmin"][1, d] = -ic_mva/data["baseMVA"]
                extradata["gen"]["4"]["cost"][d] = [-1*(da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me),0]

            #Intra ahead UK
            #source generator
                #=extradata["gen"]["7"]["pmax"][1, d] = id*ic_mva/data["baseMVA"]
                extradata["gen"]["7"]["pmin"][1, d] = 0
                extradata["gen"]["7"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]
            #load generator
                extradata["gen"]["8"]["pmax"][1, d] = -id*ic_mva/data["baseMVA"]
                extradata["gen"]["8"]["pmin"][1, d] = -id*ic_mva/data["baseMVA"]
                extradata["gen"]["8"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]=#

                extradata["gen"]["5"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["5"]["pmin"][1, d] = 0

                #extradata["gen"]["3"]["cost"][d] = [0.0]
                wnd_cost=max(0.0,da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me,da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me)
                #wnd_cost=min(0.0,21.0,1.0)
                extradata["gen"]["5"]["cost"][d] = [0.0,0.0]

            #Day ahead BE
            #source generator
                extradata["gen"]["6"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["6"]["pmin"][1, d] = 0
                extradata["gen"]["6"]["cost"][d] = [-1*(da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me),0]
            #load generator
                extradata["gen"]["7"]["pmax"][1, d] = 0
                extradata["gen"]["7"]["pmin"][1, d] = -owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["7"]["cost"][d] = [-1*(da*gencost_beuk[1,d]/e2me+id*gencostid_beuk[1,d]/e2me),0]

            #Intra ahead BE
            #source generator
                #=extradata["gen"]["12"]["pmax"][1, d] = id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["12"]["pmin"][1, d] = 0
                extradata["gen"]["12"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]
            #load generator
                extradata["gen"]["13"]["pmax"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["13"]["pmin"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                extradata["gen"]["13"]["cost"][d] = [gencostid_beuk[1,d]/e2me,0]=#


                #Day ahead UK
                #source generator
                    extradata["gen"]["8"]["pmax"][1, d] = owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["8"]["pmin"][1, d] = 0
                    extradata["gen"]["8"]["cost"][d] = [-1*(da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me),0]
                #load generator
                    extradata["gen"]["9"]["pmax"][1, d] = 0
                    extradata["gen"]["9"]["pmin"][1, d] = -owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["9"]["cost"][d] = [-1*(da*gencost_beuk[2,d]/e2me+id*gencostid_beuk[2,d]/e2me),0]

                #Intra ahead UK
                #source generator
                    #=extradata["gen"]["16"]["pmax"][1, d] = id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["16"]["pmin"][1, d] = 0
                    extradata["gen"]["16"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]
                #load generator
                    extradata["gen"]["17"]["pmax"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["17"]["pmin"][1, d] = -id*owpp_mva/data["baseMVA"]*windgenprofile_beuk[1, d]
                    extradata["gen"]["17"]["cost"][d] = [gencostid_beuk[2,d]/e2me,0]=#
    end
    return extradata
end

function create_profile_UKBE_wind_only_wstrg(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses)
    e2me=1000000/data["baseMVA"]#into ME/PU
    #e2me=1#into ME/PU
    extradata = Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    extradata["storage"] = Dict{String,Any}()
    extradata["ne_storage"] = Dict{String,Any}()
    #generator extra data
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    #existing storage
    for (b, bat) in data["storage"]
        extradata["storage"][b] = Dict{String,Any}()
        extradata["storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    #candidate storage
    for (b, bat) in data["ne_storage"]
        extradata["ne_storage"][b] = Dict{String,Any}()
        extradata["ne_storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        if (nFlow_losses[1, d]>=0)#net flow is from BE to UK
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]

            #extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["1"]["pmin"][1, d] = 0
            ##############

            #extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]]
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]
        #load generator

            extradata["gen"]["2"]["pmax"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            #extradata["gen"]["2"]["pmax"][1, d] = 0
            ##############

            #extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]#NOTE HERE
            extradata["gen"]["2"]["pmin"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            ##############

            #extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d] * -1]
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
        else#net flow is from UK to BE
            #load generator

            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            #extradata["gen"]["1"]["pmax"][1, d] = 0
            ##############

            #extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]#NOTE HERE
            extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            ##############

            #extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]
        #source generator
            extradata["gen"]["2"]["pmax"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]

            #extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["2"]["pmin"][1, d] = 0
            ##############

            #extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]]
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
        end
        extradata["gen"]["3"]["pmax"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]

        #extradata["gen"]["3"]["pmin"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        extradata["gen"]["3"]["pmin"][1, d] = 0
        ##############

        #extradata["gen"]["3"]["cost"][d] = [0.0]
        wnd_cost=min(0.0,gencost_beuk[1,d],gencost_beuk[2,d])
        #wnd_cost=min(0.0,21.0,1.0)
        extradata["gen"]["3"]["cost"][d] = [wnd_cost/e2me,0.0]

        extradata["gen"]["4"]["pmax"][1, d] = 0.0
        #NOTE changed
        ##############
        extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        #extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d] * abs(nFlow_losses[1, d])
        #extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
        extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]

        cost_abs=min(gencostid_beuk[1,d],gencostid_beuk[2,d])/e2me
        cost_inj=max(gencostid_beuk[1,d],gencostid_beuk[2,d])/e2me
        for (b, bat) in extradata["storage"]
            bat["cost_inj"][d] = 0#cost_inj
            bat["cost_abs"][d] = 0#cost_abs
        end
        #candidate storage
        for (b, bat) in extradata["ne_storage"]
            bat["cost_inj"][d] = 0#cost_inj
            bat["cost_abs"][d] = 0#cost_abs
        end
    end
    return extradata
end

#adding storage costs
function create_profile_UKBE_wstrg(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses)
    e2me=1000000/data["baseMVA"]#into ME/PU
    #e2me=1#into ME/PU
    extradata = Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    extradata["storage"] = Dict{String,Any}()
    extradata["ne_storage"] = Dict{String,Any}()
    #generator extra data
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    #existing storage
    for (b, bat) in data["storage"]
        extradata["storage"][b] = Dict{String,Any}()
        extradata["storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    #candidate storage
    for (b, bat) in data["ne_storage"]
        extradata["ne_storage"][b] = Dict{String,Any}()
        extradata["ne_storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        if (nFlow_losses[1, d]>=0)#net flow is from BE to UK
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]

            #extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["1"]["pmin"][1, d] = 0
            ##############

            #extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]]
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]
        #load generator

            extradata["gen"]["2"]["pmax"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            #extradata["gen"]["2"]["pmax"][1, d] = 0
            ##############

            #extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]#NOTE HERE
            extradata["gen"]["2"]["pmin"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            ##############

            #extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d] * -1]
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
        else#net flow is from UK to BE
            #load generator

            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            #extradata["gen"]["1"]["pmax"][1, d] = 0
            ##############

            #extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]#NOTE HERE
            extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            ##############

            #extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]
        #source generator
            extradata["gen"]["2"]["pmax"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]

            #extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["2"]["pmin"][1, d] = 0
            ##############

            #extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]]
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]
        end
        extradata["gen"]["3"]["pmax"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]

        #extradata["gen"]["3"]["pmin"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        extradata["gen"]["3"]["pmin"][1, d] = 0
        ##############

        #extradata["gen"]["3"]["cost"][d] = [0.0]
        wnd_cost=min(0.0,gencost_beuk[1,d],gencost_beuk[2,d])
        #wnd_cost=min(0.0,21.0,1.0)
        extradata["gen"]["3"]["cost"][d] = [wnd_cost/e2me,0.0]

        extradata["gen"]["4"]["pmax"][1, d] = 0.0
        #NOTE changed
        ##############
        extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        #extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d] * abs(nFlow_losses[1, d])
        #extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
        extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d]/e2me,0]

        extradata["gen"]["5"]["pmax"][1, d] = 0.0
        #NOTE changed
        ##############
        extradata["gen"]["5"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d]
        #extradata["gen"]["5"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[1, d] * abs(nFlow_losses[1, d])
        #extradata["gen"]["5"]["cost"][d] = [gencost_beuk[2,d] * -1.0]
        extradata["gen"]["5"]["cost"][d] = [gencost_beuk[2,d]/e2me,0]

        cost_abs=min(gencostid_beuk[1,d],gencostid_beuk[2,d])/e2me
        cost_inj=max(gencostid_beuk[1,d],gencostid_beuk[2,d])/e2me
        for (b, bat) in extradata["storage"]
            bat["cost_inj"][d] = 0#cost_inj
            bat["cost_abs"][d] = 0#cost_abs
        end
        #candidate storage
        for (b, bat) in extradata["ne_storage"]
            bat["cost_inj"][d] = 0#cost_inj
            bat["cost_abs"][d] = 0#cost_abs
        end
    end
    return extradata
end

#2UK
function create_profile_UKBE_2uk(number_of_hours, data, windgenprofile_beuk, gencostid_beuk, gencost_beuk, nFlow_losses)
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
        if (nFlow_losses[1, d]>=0)#net flow is from BE to UK
        #source generator
            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d]]
        #load generator
            extradata["gen"]["2"]["pmax"][1, d] = (data["gen"]["2"]["pmin"] * nFlow_losses[1, d])
            extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[2, d]#NOTE HERE
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d] * -1]
        else#net flow is from UK to BE
            #load generator
            extradata["gen"]["1"]["pmax"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]
            extradata["gen"]["1"]["pmin"][1, d] = data["gen"]["1"]["pmax"] * nFlow_losses[1, d]-data["gen"]["3"]["pmax"] * windgenprofile_beuk[2, d]#NOTE HERE
            extradata["gen"]["1"]["cost"][d] = [gencost_beuk[1,d] * -1.0]
        #source generator
            extradata["gen"]["2"]["pmax"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d]
            extradata["gen"]["2"]["pmin"][1, d] = data["gen"]["2"]["pmin"] * nFlow_losses[1, d] * 0.1
            extradata["gen"]["2"]["cost"][d] = [gencost_beuk[2,d]]
        end
        extradata["gen"]["3"]["pmax"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[2, d]
        extradata["gen"]["3"]["pmin"][1, d] = data["gen"]["3"]["pmax"] * windgenprofile_beuk[2, d]
        extradata["gen"]["3"]["cost"][d] = [0.0]

        extradata["gen"]["4"]["pmax"][1, d] = 0.0
        extradata["gen"]["4"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[2, d]
        extradata["gen"]["4"]["cost"][d] = [gencost_beuk[1,d] * -1.0]

        extradata["gen"]["5"]["pmax"][1, d] = 0.0
        extradata["gen"]["5"]["pmin"][1, d] = -1*data["gen"]["3"]["pmax"] * windgenprofile_beuk[2, d]
        extradata["gen"]["5"]["cost"][d] = [gencost_beuk[2,d] * -1.0]
    end
    return extradata
end

#Print ratings
function printEQuipValues(data)
    [println("gen"*i*" mx: "*string(d["pmax"])*" mn: "*string(d["pmin"])*" eu: "*string(d["cost"])) for (i,d) in data["gen"]]
    [println("ACbranch"*i*" mw: "*string(d["rate_a"])*" eu: "*string(d["construction_cost"])) for (i,d) in data["ne_branch"]]
    [println("DCbranch"*i*" mw: "*string(d["rateA"])*" eu: "*string(d["cost"])) for (i,d) in data["branchdc_ne"]]
    [println("Conv"*i*" mx: "*string(d["Pacmax"])*" mn: "*string(d["Pacmin"])*" eu: "*string(d["cost"])) for (i,d) in data["convdc_ne"]]
end
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

#"connection,ic_mva,ic_km,owpp_mva,owpp_km,strg,ic,ac,cnv"
function store_data2(wcfile,resultACDC, ic_mva, owpp_mva, ic_length, owpp_km)
    if (resultACDC["solution"]["nw"]["1"]["ne_branch"]["1"]["built"]>0)
        connection=-1
    elseif (resultACDC["solution"]["nw"]["1"]["ne_branch"]["2"]["built"]>0)
        connection=0
    else
        connection=1
    end
    println(wcfile,string(connection)*", "*string(resultACDC["objective"])*", "*string(resultACDC["solve_time"])*", "*string(ic_mva)*", "*string(ic_length)*", "*string(owpp_mva)*", "*string(owpp_km))
end


####################################### Clustering with bins ###################
function sample_cluster(tss_bins,tss,n)
    sample_tss=[]
    tinf=length(tss[1].time_stamp)
    loops=0
    while (length(sample_tss)<n && loops<10000)
        sample_i=rand(1:tinf)
        sample=tss[1][sample_i,:].time_stamp
        ok,tss_bins=check_sample(sample,tss,tss_bins)#working on this function
        if (ok)
            println(length(sample_tss))
            push!(sample_tss,sample)#store sample set
            loops=0
        end
        loops=loops+1
    end
    return sample_tss
end

function check_sample(sample,tss,tss_bins)
    ok=true;locations=[]
    for (key0,ts_bins) in tss_bins
        for (key1,ts_bin) in ts_bins
            if (length(filter(row -> row.time_stamp == sample, ts_bin["ts"]).time_stamp)>0)
                if (ts_bin["full"]==false)
                    push!(locations,(key0,key1))
                else
                    ok=false
                    @goto bad_sample
                end
            end
        end
    end
    if (length(locations)==length(tss_bins))
        for location in locations
            if (haskey(tss_bins[location[1]][location[2]], "chosen")==false);
                push!(tss_bins[location[1]][location[2]],"chosen"=>1);
            else
                tss_bins[location[1]][location[2]]["chosen"]=tss_bins[location[1]][location[2]]["chosen"]+1
            end
            if (tss_bins[location[1]][location[2]]["elements"]-tss_bins[location[1]][location[2]]["chosen"]==0)
                tss_bins[location[1]][location[2]]["full"]=true
            end
        end
    else
        ok=false
    end
    @label bad_sample
    return ok,tss_bins
end

function cluster_probs(tss_bins,tss,n)
    tinf=length(tss[1].time_stamp)
    for (k0,ts_bins) in tss_bins
        for (k1,ts_bin) in ts_bins
            elements=length(ts_bin["ts"].time_stamp)==0 ? 0 : ceil(Int64,n*length(ts_bin["ts"].time_stamp)/tinf)
            push!(ts_bin,"elements"=>elements)
            if (elements>0); push!(ts_bin,"full"=>false);else push!(ts_bin,"full"=>true);end
        end
    end
    return tss_bins
end
#tss is an array of time series arays, n is the desired sample size
function cluster_ts(tss,n)
    number_tss=length(tss)
    tss_bins=Dict{String,Any}()
    for (key,ts) in enumerate(tss)
        push!(tss_bins,string(key)=>ts_binify(ts))#sort into bins
    end
    tss_bins=cluster_probs(tss_bins,tss,n)
    return tss_bins
end

function ts_binify(ts)
    ts_bins=Dict{String,Any}()
    nrm_ts=normalize_ts(deepcopy(ts))
    nm_bins=10
    for b=1/nm_bins:1/nm_bins:1
        for p in eachrow(nrm_ts)
            if (p.price<=b && p.price>(b-1/nm_bins))
                if (haskey(ts_bins, string(b))==false);push!(ts_bins,string(b)=>Dict{String,Any}());push!(ts_bins[string(b)],"ts"=>DataFrame("time_stamp"=>[]));end
                push!(ts_bins[string(b)]["ts"],p[1:1])
            end
        end
    end
    return ts_bins
end

function normalize_ts(nrm_ts)
    mn=minimum(nrm_ts.price)
    nrm_ts.price=nrm_ts.price.-mn#set minimum to zero
    mx=maximum(nrm_ts.price)
    nrm_ts.price=nrm_ts.price./mx#normalize
    return nrm_ts
end
