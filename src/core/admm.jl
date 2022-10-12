#=Main logic for ADMM
function admm_4_AjAwAgAuAo_main(mn_data, gurobi, s)
    results_set=[]
    eps=s["eps"]#set max value of residual for convergence
    residual=Inf#Initialize residual
    #agents=["Ag","Au","Ao","Aw","Aj"]#Initilize agents
    agents=["Aall","Ao"]#Initilize agents
    push!(s,"fixed_variables" => Dict{String,Any}())#reserve memory for varible updates
    push!(s,"agent" => "")#create agent dictionary entry
    s["fixed_variables"] = admm_4_AjAwAgAuAo_intialize(mn_data["nw"], s["fixed_variables"], s["wfz"], s["genz"])#set initial values for fixed variables
    while (residual>eps)
    #for i=1:10
        results_set=[]#reserve memory for results
        for a in agents
            s["agent"]=a#set agent
            println(a)
            result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting=s)#solve sub problem
            if !(isnan(result_mip["objective"]))#check solution exists
                s["fixed_variables"] = admm_4_AjAwAgAuAo_update(result_mip["solution"]["nw"], s["fixed_variables"], a, s["wfz"])#update agent fixed variables
                push!(results_set,(a,result_mip))#store results
            else
                println("skyap!!!")#if no solution re-iterate
            end
        end
        [println(first(r)*" "*string(last(r)["objective"])) for r in results_set]#print current objective
        s["fixed_variables"], residual = dual_variable_update(s["fixed_variables"], s["beta"])#update the dual variable
        println("Residual: "*string(residual))#print residual
    #end
    end
    return last(results_set[length(results_set)])
end

################################################# Initilize agents ###############################################
function admm_4_AjAwAgAuAo_intialize(mn_data_nw, fixed_variables, wfz, genz)
    for (key,nw) in mn_data_nw
        push!(fixed_variables,key=>Dict{String,Any}())
        push!(fixed_variables[key],"baseMVA" => nw["baseMVA"])
        #initialize all gen values in dictionary to zero
        push!(fixed_variables[key],"gen" => Dict{String, Any}())
        for (key_g,g) in nw["gen"]
            int_g=parse(Int64,key_g)
            push!(fixed_variables[key]["gen"],key_g => Dict{String, Any}())
            push!(fixed_variables[key]["gen"][key_g],"qg" => 0.00)
            if (issubset(int_g, first.(genz)) && int_g < minimum(first.(wfz)))#conventional gens
                push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            elseif (issubset(int_g, first.(genz)) && int_g > maximum(first.(wfz)))#consumers
                push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            elseif (issubset(parse(Int64,key_g), first.(wfz)))#wind farms
                push!(fixed_variables[key]["gen"][key_g],"wf_pacmax" => 0.00)
                push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            end
        end
        #initialize all bus values in dictionary to zero
        push!(fixed_variables[key],"bus" => Dict{String, Any}())
        for (key_b,b) in nw["bus"]
            push!(fixed_variables[key]["bus"],key_b => Dict{String, Any}())
            push!(fixed_variables[key]["bus"][key_b],"va" => 0.00)
            push!(fixed_variables[key]["bus"][key_b],"lambda" => 1.0)
            push!(fixed_variables[key]["bus"][key_b],"imbalance" => 0.00)
        end
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"branchdc_ne" => Dict{String, Any}())
        for (key_br,br) in nw["branchdc_ne"]
            push!(fixed_variables[key]["branchdc_ne"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["branchdc_ne"][key_br],"isbuilt" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"pt" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"pf" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"fbusdc" => br["fbusdc"])
            push!(fixed_variables[key]["branchdc_ne"][key_br],"tbusdc" => br["tbusdc"])
            push!(fixed_variables[key]["branchdc_ne"][key_br],"rateA" => br["rateA"])
        end

        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"ne_branch" => Dict{String, Any}())
        for (key_br,br) in nw["ne_branch"]
            push!(fixed_variables[key]["ne_branch"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["ne_branch"][key_br],"built" => 0.00)
            push!(fixed_variables[key]["ne_branch"][key_br],"p_ne_to" => 0.00)
            push!(fixed_variables[key]["ne_branch"][key_br],"p_ne_fr" => 0.00)
            push!(fixed_variables[key]["ne_branch"][key_br],"f_bus" => br["f_bus"])
            push!(fixed_variables[key]["ne_branch"][key_br],"t_bus" => br["t_bus"])
            push!(fixed_variables[key]["ne_branch"][key_br],"rate_a" => br["rate_a"])
        end

        push!(fixed_variables[key],"branchdc" => Dict{String, Any}())
        for (key_br,br) in nw["branchdc"]
            push!(fixed_variables[key]["branchdc"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["branchdc"][key_br],"p_rateA" => 0.00)
            push!(fixed_variables[key]["branchdc"][key_br],"pt" => 0.00)
            push!(fixed_variables[key]["branchdc"][key_br],"pf" => 0.00)
            push!(fixed_variables[key]["branchdc"][key_br],"fbusdc" => br["fbusdc"])
            push!(fixed_variables[key]["branchdc"][key_br],"tbusdc" => br["tbusdc"])
        end

        push!(fixed_variables[key],"branch" => Dict{String, Any}())
        for (key_br,br) in nw["branch"]
            push!(fixed_variables[key]["branch"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["branch"][key_br],"p_rateAC" => 0.00)
            push!(fixed_variables[key]["branch"][key_br],"pt" => 0.00)
            push!(fixed_variables[key]["branch"][key_br],"pf" => 0.00)
            push!(fixed_variables[key]["branch"][key_br],"f_bus" => br["f_bus"])
            push!(fixed_variables[key]["branch"][key_br],"t_bus" => br["t_bus"])
        end
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"storage" => Dict{String, Any}())
        for (key_s,s) in nw["storage"]
            push!(fixed_variables[key]["storage"],key_s => Dict{String, Any}())
            push!(fixed_variables[key]["storage"][key_s],"qs" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"e_absmax" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"e_abs" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"qsc" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"se" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"sc" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"sd" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"ps" => 0.00)
        end
        #if needed add convdc_ne here
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"convdc" => Dict{String, Any}())
        for (key_c,c) in nw["convdc"]
            push!(fixed_variables[key]["convdc"],key_c => Dict{String, Any}())
            push!(fixed_variables[key]["convdc"][key_c],"vaf" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_pr_fr" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_dc" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"p_pacmax" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"vac" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_ac" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_tf_to" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_tf_fr" => 0.00)
        end
    end
    return fixed_variables
end

################################################# Update Agents ##################################################
function admm_4_AjAwAgAuAo_update(rez, fixed_variables, agent, wfz)
    for (n,nw) in rez
        if (agent=="Ao")#Only convex cables are used
            fixed_variables=update_branchdc(n,nw,fixed_variables)
            fixed_variables=update_branch(n,nw,fixed_variables)
            fixed_variables=update_dcconv(n,nw,fixed_variables)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Au" || agent=="Ag" || agent=="Aw")#consumpttion, conventional, wind
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Aj")#storage
            fixed_variables=update_storage(n,nw,fixed_variables)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Aall")#all agents except TSO
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            fixed_variables=update_storage(n,nw,fixed_variables)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        end
    end
    return fixed_variables
end

#update all convdc values in dictionary
function update_dcconv(key,nw,fixed_variables)
    for (key_c,c) in nw["convdc"]
        fixed_variables[key]["convdc"][key_c]["vaf"] = nw["convdc"][key_c]["vafilt"]
        fixed_variables[key]["convdc"][key_c]["pconv_pr_fr"] = nw["convdc"][key_c]["ppr_fr"]
        fixed_variables[key]["convdc"][key_c]["pconv_dc"] = nw["convdc"][key_c]["pdc"]
        fixed_variables[key]["convdc"][key_c]["p_pacmax"] = nw["convdc"][key_c]["p_pacmax"]
        fixed_variables[key]["convdc"][key_c]["vac"] = nw["convdc"][key_c]["vaconv"]
        fixed_variables[key]["convdc"][key_c]["pconv_ac"] = nw["convdc"][key_c]["pconv"]
        fixed_variables[key]["convdc"][key_c]["pconv_tf_to"] = nw["convdc"][key_c]["ptf_to"]
        fixed_variables[key]["convdc"][key_c]["pconv_tf_fr"] = -1*nw["convdc"][key_c]["ptf_to"]#nw["convdc"][key_c]["pgrid"]
    end
    return fixed_variables
end

#updates conventional and wind generation
function update_genz(key,nw,fixed_variables,wfz)
    for (key_g,g) in nw["gen"]
        fixed_variables[key]["gen"][key_g]["qg"] = nw["gen"][key_g]["qg"]
        fixed_variables[key]["gen"][key_g]["pg"] = nw["gen"][key_g]["pg"]
        if (issubset(parse(Int64,key_g), first.(wfz)))
            fixed_variables[key]["gen"][key_g]["wf_pacmax"] = nw["gen"][key_g]["wf_pacmax"]
        end
    end
    return fixed_variables
end

#updates recorded imbalance for each node
function update_imbalance(key,nw,fixed_variables)
    #update all storage values in dictionary
    for (key_b,b) in nw["bus"]
        fixed_variables[key]["bus"][key_b]["imbalance"] = nw["bus"][key_b]["imbalance"]
    end
    return fixed_variables
end

######### Cables
#update convex cable dictionary
function update_branchdc(key,nw,fixed_variables)
    for (key_br,br) in nw["branchdc"]
        fixed_variables[key]["branchdc"][key_br]["p_rateA"] = nw["branchdc"][key_br]["p_rateA"]
        fixed_variables[key]["branchdc"][key_br]["pt"] = nw["branchdc"][key_br]["pt"]
        fixed_variables[key]["branchdc"][key_br]["pf"] = -1*nw["branchdc"][key_br]["pt"]
    end
    return fixed_variables
end

function update_branch(key,nw,fixed_variables)
    for (key_br,br) in nw["branch"]
        fixed_variables[key]["branch"][key_br]["p_rateAC"] = nw["branch"][key_br]["p_rateAC"]
        fixed_variables[key]["branch"][key_br]["pt"] = nw["branch"][key_br]["pt"]
        fixed_variables[key]["branch"][key_br]["pf"] = -1*nw["branch"][key_br]["pt"]
    end
    return fixed_variables
end

#update all storage values in dictionary
function update_storage(key,nw,fixed_variables)
    for (key_s,s) in nw["storage"]
        fixed_variables[key]["storage"][key_s]["qs"] = nw["storage"][key_s]["qs"]
        fixed_variables[key]["storage"][key_s]["e_absmax"] = nw["storage"][key_s]["e_absmax"]
        fixed_variables[key]["storage"][key_s]["e_abs"] = nw["storage"][key_s]["e_abs"]
        fixed_variables[key]["storage"][key_s]["qsc"] = nw["storage"][key_s]["qsc"]
        fixed_variables[key]["storage"][key_s]["se"] = nw["storage"][key_s]["se"]
        fixed_variables[key]["storage"][key_s]["sc"] = nw["storage"][key_s]["sc"]
        fixed_variables[key]["storage"][key_s]["sd"] = nw["storage"][key_s]["sd"]
        fixed_variables[key]["storage"][key_s]["ps"] = nw["storage"][key_s]["ps"]
    end
    return fixed_variables
end
################################################# Dual update ####################################################
function dual_variable_update(fixed_variables, beta)
    residuals=[]
    for (n,nw) in fixed_variables
        for (i_b,b) in nw["bus"]
            lambda0=deepcopy(b["lambda"])
            b["lambda"]=lambda0+beta*b["imbalance"]
            push!(residuals,abs(b["lambda"]-lambda0))
        end
    end
    return fixed_variables, maximum(residuals)
end

################################################# fixing agent values ############################################
function fix_variables(pm)
    for n in _PM.nw_ids(pm)
        if (pm.setting["agent"]=="Ao")
            fix_storage(pm, n)#set storage vars
            fix_genz(pm, n)#set gen vars
            fix_wind(pm, n)#set wind vars
            fix_loadz(pm, n)#set loads
        elseif (pm.setting["agent"]=="Aw")
            fix_storage(pm, n)#set storage vars
            fix_genz(pm, n)#set gen vars
            fix_loadz(pm, n)#set loads
            fix_branchesdc(pm, n)#set branch vars
            fix_branchesac(pm, n)
            fix_convdc(pm,n)#set convdc vars
        elseif (pm.setting["agent"]=="Aj")
            fix_wind(pm, n)#set wind vars
            fix_genz(pm, n)#set gen vars
            fix_loadz(pm, n)#set loads
            fix_branchesdc(pm, n)#set branch vars
            fix_branchesac(pm, n)
            fix_convdc(pm,n)#set convdc vars
        elseif (pm.setting["agent"]=="Ag")
            fix_storage(pm, n)#set storage vars
            fix_wind(pm, n)#set wind vars
            fix_loadz(pm, n)#set loads
            fix_branchesdc(pm, n)#set branch vars
            fix_branchesac(pm, n)
            fix_convdc(pm,n)#set convdc vars
        elseif (pm.setting["agent"]=="Au")
            fix_storage(pm, n)#set storage vars
            fix_wind(pm, n)#set wind vars
            fix_genz(pm, n)#set gen vars
            fix_branchesdc(pm, n)#set branch vars
            fix_branchesac(pm, n)
            fix_convdc(pm,n)#set convdc vars
        elseif (pm.setting["agent"]=="Aall")
            fix_branchesdc(pm, n)#set branch vars
            fix_branchesac(pm, n)
            fix_convdc(pm,n)#set convdc vars
        end
    end
end


function fix_storage(pm, n)#set storage vars
    for (i,s) in pm.setting["fixed_variables"][string(n)]["storage"]
        for (kv,v) in s
            storage = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
            JuMP.fix(storage,v,force=true)
        end
    end
end

function fix_genz(pm, n)#set conv genz

    for (i,g) in pm.setting["fixed_variables"][string(n)]["gen"]
        i_int=parse(Int64,i)
        if (i_int<minimum(first.(pm.setting["wfz"])))
            #println("genz: "*i)
            for (kv,v) in g
                gen = _PM.var(pm, n, Symbol(kv), i_int)
                JuMP.fix(gen,v,force=true)
            end
        end
    end
end

function fix_loadz(pm, n)#set loads

    for (i,g) in pm.setting["fixed_variables"][string(n)]["gen"]
        i_int=parse(Int64,i)
        if (i_int>maximum(first.(pm.setting["wfz"])))
            #println("loadz: "*i)
            for (kv,v) in g
                #println(kv*" value "*string(v))
                gen = _PM.var(pm, n, Symbol(kv), i_int)
                JuMP.fix(gen,v,force=true)
            end
        end
    end
end

function fix_wind(pm, n)#set wind

    for (i,g) in pm.setting["fixed_variables"][string(n)]["gen"]
        i_int=parse(Int64,i)
        if (issubset(i_int, first.(pm.setting["wfz"])))
            for (kv,v) in g
                gen = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
                JuMP.fix(gen,v,force=true)
            end
        end
    end
end

function fix_convdc(pm,n)
    #set convdc
    for (i,cv) in pm.setting["fixed_variables"][string(n)]["convdc"]
        for (kv,v) in cv
            #if !(issubset([kv],[]))
                conv = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
                JuMP.fix(conv,v,force=true)
            #end
        end
    end
end


function fix_branchesdc(pm, n)#set branch vars
    #set convdc
    if haskey(_PM.var(pm, n), :p_dcgrid)
        p_dcgrid = _PM.var(pm, n, :p_dcgrid)
        if !isempty(p_dcgrid)
            brchs=[]
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_from)
                JuMP.fix(p_dcgrid[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branchdc"][string(l)]["pf"],force=true)
            end
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_to)
                JuMP.fix(p_dcgrid[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branchdc"][string(l)]["pt"],force=true)
            end
        end
    end
end

function fix_branchesac(pm, n)#set branch vars
    #set convdc
    if haskey(_PM.var(pm, n), :p)
        _p = _PM.var(pm, n, :p)
        if !isempty(_p)
            brchs=[]
            for (l,i,j) in _PM.ref(pm, n, :arcs_from)
                JuMP.fix(_p[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branch"][string(l)]["pf"],force=true)
            end
            for (l,i,j) in _PM.ref(pm, n, :arcs_to)
                JuMP.fix(_p[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branch"][string(l)]["pt"],force=true)
            end
        end
    end
end
=#