################################################ start cable crossing  constraints pre-rocessing ################################# 

function all_cable_intersections(s,data)
    intersections=allIntersects(s, data)
    ac_cable_intersections=ac_intersections(intersections["AC"], data["ne_branch"])
    dc_cable_intersections=dc_intersections(intersections["DC"], data["branchdc_ne"], data["convdc"])
    acdc_cable_intersections=acdc_intersections(intersections["ACDC"], data["ne_branch"], data["branchdc_ne"], data["convdc"])
    return Dict("ac"=>ac_cable_intersections, "dc"=>dc_cable_intersections, "acdc"=>acdc_cable_intersections)  
end  

function dc_intersections(dc_intersections, dc_cables, convdc) 
    #if DC cables are within a collection circuit - this will be required
    dc_cable_intersections=DataFrames.DataFrame()
    return dc_cable_intersections
end


function acdc_intersections(acdc_intersections, ac_cables, dc_cables, convdc) 
    #if DC cables are within a collection circuit - this will be required
    acdc_cable_intersections=DataFrames.DataFrame()
    return acdc_cable_intersections
end

function ac_intersections(ac_intersections, ac_cables) 
    ac_cable_intersections=DataFrames.DataFrame()
    cables_1=[];cables_2=[]
    for _row in eachrow(ac_intersections)
        cable_1=[];cable_2=[]
        for (key_ac,ac_cable) in ac_cables 
            if (ac_cable["f_bus"]==_row[:fr_ac1] && ac_cable["t_bus"]==_row[:to_ac1])
                push!(cable_1,key_ac);end
            if (ac_cable["f_bus"]==_row[:fr_ac2] && ac_cable["t_bus"]==_row[:to_ac2])
                push!(cable_2,key_ac);end
        end
        unique!(cable_1)
        unique!(cable_2)
        push!(cables_1,cable_1)
        push!(cables_2,cable_2)
    end
    ac_cable_intersections[!,:ac_cables_1]=cables_1
    ac_cable_intersections[!,:ac_cables_2]=cables_2
    return ac_cable_intersections
end

#Given three collinear points p, q, r, the function checks if point q lies on line segment 'pr'
function onSegment(p, q, r)
        if (first(q) <= max(first(p), first(r)) && first(q) >= min(first(p), first(r)) && last(q) <= max(last(p), last(r)) && last(q) >= min(last(p), last(r)))
           return true;end
      
        return false;
end
      

function orientation(p, q, r)
        #See https://www.geeksforgeeks.org/orientation-3-ordered-points/ for details of below formula.
        val = (last(q) - last(p)) * (first(r) - first(q)) - (first(q) - first(p)) * (last(r) - last(q));
      
        if (val == 0) return 0;end  # collinear
      
        return val>0 ? 1 : 2; # clock or counterclock wise
end
    

function doIntersect(p1, q1, p2, q2)
    #If a terminal node - matches there is no overlap
    if (p1==p2 || p1==q2 || q1==p2 || q1==q2)
        return false;end
    
        #Find the four orientations needed for general and special cases
    o1 = orientation(p1, q1, p2);
    o2 = orientation(p1, q1, q2);
    o3 = orientation(p2, q2, p1);
    o4 = orientation(p2, q2, q1);
    
    #General case
    if (o1 != o2 && o3 != o4)
        return true;end
    
    #Special Cases p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;end
    
    #p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;end
    
    #p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;end
    
    #p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;end
    
    return false; # Doesn't fall in any of the above cases
end

function allIntersects(s,data)
    df_map=map_Of_Connections_ACDCNTC(data, s)
    intersections=Dict("DC"=>DataFrames.DataFrame(:fr_dc1=>[],:to_dc1=>[],:fr_dc2=>[],:to_dc2=>[]),"AC"=>DataFrames.DataFrame(:fr_ac1=>[],:to_ac1=>[],:fr_ac2=>[],:to_ac2=>[]),"ACDC"=>DataFrames.DataFrame(:fr_ac1=>[],:to_ac1=>[],:fr_dc2=>[],:to_dc2=>[]))
    for (i1,_row1) in enumerate(eachrow(df_map))
        for (i2,_row2) in enumerate(eachrow(df_map))
            if (i1<i2)
                acdc1=string(_row1[:type])
                fr1=(_row1[:fr_node])
                to1=(_row1[:to_node])
                p1=(_row1[:lat_fr],_row1[:long_fr])
                q1=(_row1[:lat_to],_row1[:long_to])

                acdc2=string(_row2[:type])
                fr2=(_row2[:fr_node])
                to2=(_row2[:to_node])
                p2=(_row2[:lat_fr],_row2[:long_fr])
                q2=(_row2[:lat_to],_row2[:long_to])

                if (acdc1!="NTC" && acdc2!="NTC")
                    test_result=doIntersect(p1, q1, p2, q2)
                    if (test_result==true)
                        if (acdc1=="AC" && acdc2=="AC")
                            push!(intersections["AC"],[fr1,to1,fr2,to2])
                        elseif (acdc1=="DC" && acdc2=="DC")
                            push!(intersections["DC"],[fr1,to1,fr2,to2])
                        elseif (acdc1=="AC")
                            push!(intersections["ACDC"],[fr1,to1,fr2,to2])
                        else
                            push!(intersections["ACDC"],[fr2,to2,fr1,to1])
                        end
                    else
                    end
                end
            end
        end
    end
    return intersections
end
################################################ end cable crossing  constraints pre-rocessing ################################# 

################################################ max strings per oss constraints pre-rocessing ################################# 
function oss_connections(s,data)
    ac_connections=data["ne_branch"];
    dc_connections=data["branchdc_ne"];
    cables_connected_to_osss=Dict{Int64,Any}()
    
    for oss_node in s["oss_nodes"] 
        
        #find all ac cables connected to the oss
        push!(cables_connected_to_osss,oss_node=>Dict("ac"=>[],"dc"=>[]))
        for (ac_key,ac_connection) in ac_connections
            if (ac_connection["f_bus"]==oss_node || ac_connection["t_bus"]==oss_node)
                push!(cables_connected_to_osss[oss_node]["ac"],ac_key)
            end
        end

        #find the dc node connected to the ac node
        busdc_i=0
        for (conv_key, convdc) in data["convdc"]
            if (convdc["busac_i"]==oss_node)
                busdc_i=convdc["busdc_i"]
            end
        end

        #find all dc cables connected to the oss
        for (dc_key,dc_connection) in dc_connections        
            if (dc_connection["fbusdc"]==busdc_i || dc_connection["tbusdc"]==busdc_i)
                push!(cables_connected_to_osss[oss_node]["dc"],dc_key)
            end
        end
    end

    return cables_connected_to_osss
end

function oss_string_feeders(s,data)
    onshore_nodes=filter(:type=>x->x!=0,s["nodes"])[!,:node]
    ac_connections=data["ne_branch"];
    dc_connections=data["branchdc_ne"];
    cables_connected_to_osss=Dict{Int64,Any}()
    
    for oss_node in s["oss_nodes"] 
        
        #find all ac cables connected to the oss
        push!(cables_connected_to_osss,oss_node=>Dict("ac"=>[],"dc"=>[]))
        for (ac_key,ac_connection) in ac_connections
            if ((ac_connection["f_bus"]==oss_node || ac_connection["t_bus"]==oss_node) && !issubset([ac_connection["f_bus"]],onshore_nodes) && !issubset([ac_connection["t_bus"]],onshore_nodes))
                push!(cables_connected_to_osss[oss_node]["ac"],ac_key)
            end
        end

        #find the dc node connected to the ac node
        busdc_i=0
        for (conv_key, convdc) in data["convdc"]
            if (convdc["busac_i"]==oss_node)
                busdc_i=convdc["busdc_i"]
            end
        end

        #find all dc cables connected to the oss
        for (dc_key,dc_connection) in dc_connections        
            if ((dc_connection["fbusdc"]==busdc_i || dc_connection["tbusdc"]==busdc_i) && !issubset([dc_connection["fbusdc"]],onshore_nodes) && !issubset([dc_connection["tbusdc"]],onshore_nodes))
                push!(cables_connected_to_osss[oss_node]["dc"],dc_key)
            end
        end
    end

    return cables_connected_to_osss
end
################################################ end max strings per oss constraints pre-rocessing ################################# 

################################################ max connections per turbine constraints pre-rocessing ################################# 
function turbine_connections(s,data)
    ac_connections=data["ne_branch"];
    dc_connections=data["branchdc_ne"];
    cables_connected_to_turbines=Dict{Int64,Any}()

    nodes=filter(:type=>x->x==0,s["nodes"])
    filter!(:node=>x->!issubset([x],s["oss_nodes"]),nodes)
    push!(s,"turbine_nodes"=>nodes[!,:node])
    for turbine_node in nodes[!,:node] 
        
        #find all ac cables connected to the oss
        push!(cables_connected_to_turbines,turbine_node=>Dict("ac"=>[],"dc"=>[]))
        for (ac_key,ac_connection) in ac_connections
            if (ac_connection["f_bus"]==turbine_node || ac_connection["t_bus"]==turbine_node)
                push!(cables_connected_to_turbines[turbine_node]["ac"],ac_key)
            end
        end

        #find the dc node connected to the ac node
        busdc_i=0
        for (conv_key, convdc) in data["convdc"]
            if (convdc["busac_i"]==turbine_node)
                busdc_i=convdc["busdc_i"]
            end
        end

        #find all dc cables connected to the oss
        for (dc_key,dc_connection) in dc_connections        
            if (dc_connection["fbusdc"]==busdc_i || dc_connection["tbusdc"]==busdc_i)
                push!(cables_connected_to_turbines[turbine_node]["dc"],dc_key)
            end
        end
    end
    return cables_connected_to_turbines
end

#to avoid loops only one input cable to turbine is allowed
function turbine_connections_input(s,data)
    ac_connections=data["ne_branch"];
    dc_connections=data["branchdc_ne"];
    cables_connected_to_turbines=Dict{Int64,Any}()

    nodes=filter(:type=>x->x==0,s["nodes"])
    filter!(:node=>x->!issubset([x],s["oss_nodes"]),nodes)
    push!(s,"turbine_nodes"=>nodes[!,:node])
    for turbine_node in nodes[!,:node] 
        
        #find all ac cables connected to the oss
        push!(cables_connected_to_turbines,turbine_node=>Dict("ac"=>[],"dc"=>[]))
        for (ac_key,ac_connection) in ac_connections
            if (ac_connection["t_bus"]==turbine_node)
                push!(cables_connected_to_turbines[turbine_node]["ac"],ac_key)
            end
        end

        #find the dc node connected to the ac node
        busdc_i=0
        for (conv_key, convdc) in data["convdc"]
            if (convdc["busac_i"]==turbine_node)
                busdc_i=convdc["busdc_i"]
            end
        end

        #find all dc cables connected to the oss
        for (dc_key,dc_connection) in dc_connections        
            if (dc_connection["tbusdc"]==busdc_i)
                push!(cables_connected_to_turbines[turbine_node]["dc"],dc_key)
            end
        end
    end
    return cables_connected_to_turbines
end