gen_gps_location=(51.71997,2.74)
pcc_gps_location=(51.32694, 3.18361)
gen_locations=[(482039.71495261084, 5729925.21)]
pcc_locations=[(94894.12267,5702260.96)]
x0=482039.71495261084
x1=94894.12267
y0=5729925.21
y1=5702260.96
lof_rotateAxis(482039.71495261084,94894.12267,5729925.21,5702260.96,gen_locations,pcc_locations)
xy=lof_rotatePnt(x1,y1,theta)
ocean.offset

lof_cartesian2gps(ocean.gens,base)
ocean.gens[4].coord.x=18.215599597767707
ocean.gens[4].coord.y=0
function lof_layoutEez(cnt)

    ocean=eez()#build the eez in the ocean
    ocean.pccs=lof_layPccs()#add the gps of PCCs
    println("PCCs positioned at:")
    for value in ocean.pccs
        print(value.num)
        print(" - ")
        println(value.gps)
    end

    ocean.gens=lof_layGens(ocean)#add the gps of OWPPs
    println("OWPPs positioned at:")
    for value in ocean.gens
        print(value.num)
        print(" - ")
        println(value.gps)
    end
    base=lof_bseCrd(ocean)#find base coordinates
    print("Base coordinates are: ")
    println(base)
    lof_gps2cartesian(ocean.gens,base)#projects owpps onto cartesian plane
    lof_gps2cartesian(ocean.pccs,base)#projects pccs onto cartesian plane
    println("GPS coordinates projected onto cartesian plane.")
    lof_transformAxis(ocean)
    println("Axis transformed.")
    lof_osss(ocean,cnt)#add all osss within boundary
    println("OSSs positioned at:")
    for value in ocean.osss
        print(value.num)
        print(" - ")
        println(value.coord)
    end
    lof_layoutEez_arcs(ocean,cnt)
    #ppf_printOcnXY(ocean)
    return ocean
end

function lof_gensOrder(gens,ocn,num)
    lnths=Array{Float64,1}()
    ordrdGens=Array{node,1}()
    for gen in gens
        pcc_close=lof_xClosestPcc(gen,ocn.pccs)
        push!(lnths,lof_pnt2pnt_dist(gen.coord,pcc_close.coord))
    end
    for lps=1:1:length(lnths)
        push!(ordrdGens,gens[findmin(lnths)[2]])
        ordrdGens[length(ordrdGens)].num=deepcopy(num+lps)
        lnths[findmin(lnths)[2]]=Inf
    end
    gens = ordrdGens
    return gens
end

#returns the hypotenuse distance between 2 cartesian points
#minimum distance for a path is 1km
function lof_pnt2pnt_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
    if hyp < 1
        hyp=1
        #println("Arc distance is less than 1km, set to 1km.")
    end
    return hyp
end

#find closest PCC
function lof_xClosestPcc(i,pc)
    x_close=node()
    lnths=Array{Float64,1}()
    for j=1:length(pc)
        push!(lnths,lof_pnt2pnt_dist(i.coord,pc[j].coord))
    end
    mn=findmin(lnths)
    x_close=pc[mn[2]]
    return x_close
end

function lod_pccGps()
    pcc=Array{Tuple,1}()
    ##### Belgium #####
        push!(pcc,(2.939692,51.239737))
        push!(pcc,(3.183611,51.32694))
    return pcc
end
function lod_pccKv()
    return 220.0
end

#set offshore transmission voltage
function lod_ossKv()
    return 220.0
end

#set collector voltage
function lod_cncsKv()
    return 66.0
end
##############################
########## Gens ##############

function lod_gensGps()
    c=Array{Tuple,1}()
    wnd=Array{String,1}()
    p=Array{Float64,1}()
##### Belgium #####
    
        #Norther
        push!(c,(3.015833,51.52806))
        push!(p,250.0)
        push!(wnd,"Norther")
        #Thornton
        push!(c,((2.97+2.919972)/2,(51.56+51.53997)/2))
        push!(p,250.0)
        push!(wnd,"Thornton")
        #Rentel
        push!(c,(2.939972,51.59))
        push!(p,250.0)
        push!(wnd,"Rentel")
        #Northwind
        push!(c,(2.900972,51.61897))
        push!(p,250.0)
        push!(wnd,"Northwind")
        #Seastar
        push!(c,(2.859972,51.63))
        push!(p,250.0)
        push!(wnd,"Seastar")
        #Nobelwind/Belwind
        push!(c,((2.819972+2.799972)/2,(51.664+51.67)/2))
        push!(p,250.0)
        push!(wnd,"Nobelwind")
        #Northwester
        push!(c,(2.757,51.68597))
        push!(p,250.0)
        push!(wnd,"Northwester")
        #Mermaid
        push!(c,(2.74,51.71997))
        push!(p,250.0)
        push!(wnd,"Mermaid")
    return c,p,wnd
end

###################################################################
################################################################################
###################### mapping args ############################################
################################################################################
#struct used to pass mapping connection arguments
mutable struct control
    xrad::Bool
    neib1::Bool
    neib3::Bool
    xradPcc::Bool
    xradHlf::Bool
    xXrad::Array{Int64,1}
    xXneib1::Array{Int64,1}
    xXneib3::Array{Int64,1}
end
control()=control(false,false,false,false,false,[],[],[])
###################################################################
mutable struct xy
    x::Float64
    y::Float64
end
xy()=xy(69.69,69.69)
###################################################################
mutable struct gps
    lat::Float64
    lng::Float64
end
gps()=gps(69.69,69.69)
################################################################################
###################### Nodes ###################################################
################################################################################
mutable struct node
    gps::gps
    name::String
    coord::xy
    mva::Float64
    kv::Float64
    mvas::Array{Float64}
    wnds::Array{String}
    num::Int64
    id::String
    upstrm::Int64
    dwnstrm::Int64
end
node()=node(gps(),"colruyt",xy(),69.69,69.69,[],[],69,"sixty-nine",69,69)
####################################################################
##################### Arcs #########################################
####################################################################
mutable struct arc
    head::node
    tail::node
    lngth::Float64
    mva::Float64#Only used for disaplaying solution
end
arc()=arc(node(),node(),69.69,69.69)
###################################################################
mutable struct eez
    osss::Array{node}
    gens::Array{node}
    pccs::Array{node}
    gOarcs::Array{arc}
    oOarcs::Array{arc}
    oParcs::Array{arc}
    gParcs::Array{arc}
    angle::Float64
    offset::Float64
    mnGap::Float64
    oOcbls::Array{Tuple}
    oPcbls::Array{Tuple}
    oPXcbls::Array{Tuple}
    gOcbls::Array{Tuple}
    gPcbls::Array{Tuple}
    dcCbls::Array{Tuple}
end
eez()=eez([],[],[],[],[],[],[],69.69,69.69,69.69,[],[],[],[],[],[])
################################################################################
############################ GPS to cartesian transform ########################
################################################################################
#sets the gps coords that used as reference coords
function lof_bseCrd(ocean)
    base=gps()
    #for india (type layouts)
    if ocean.gens[length(ocean.gens)].gps.lat < ocean.pccs[length(ocean.pccs)].gps.lat
        base.lat=ocean.gens[length(ocean.gens)].gps.lat#base lat
        base.lng=ocean.gens[length(ocean.gens)].gps.lng#base long
    #for belgium (type layouts)
    elseif ocean.pccs[length(ocean.pccs)].gps.lat < ocean.gens[length(ocean.gens)].gps.lat
        base.lat=ocean.pccs[length(ocean.pccs)].gps.lat#base lat
        base.lng=ocean.pccs[length(ocean.pccs)].gps.lng#base long
    else
        error("No proper base coordinates system established!")
    end
    return base
end

#calculates lengths based on latitude
#as lattitude changes number of km should be updated
function lof_gps2cartesian(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.coord.x=lof_deg2lgth(value.gps.lng-base.lng,lof_lg1deg(value.gps.lat,lnthLT))
        value.coord.y=lof_deg2lgth(value.gps.lat-base.lat,lnthLT)
    end
end
ocn=ocean
#rotates and slides cartesian axis
function lof_transformAxis(ocn)
    offset=lof_rotateAxis(ocn)
    lof_slideAxis(ocn,offset)
    num=length(ocn.pccs)
    ocn.gens=lof_gensOrder(ocn.gens,ocn,num)
    ocn.mnGap=lof_mnGap(ocn.gens)
end

#finds angle to rotate and applies to owpps and pccs
#rotates axis to align n-s with y
function lof_rotateAxis(ocn)
    theta=atan((ocn.pccs[length(ocn.pccs)].coord.x-ocn.gens[length(ocn.gens)].coord.x)/(ocn.gens[length(ocn.gens)].coord.y-ocn.pccs[length(ocn.pccs)].coord.y))
    ocn.angle=theta
    offset=0.0
    offset=lof_rotateGroup(ocn.gens,theta,offset)
    offset=lof_rotateGroup(ocn.pccs,theta,offset)
    ocn.offset=offset
    return offset
end

#loops through to apply rotations for a specified group
function lof_rotateGroup(locations,theta,os)
    for value in locations
        xy=lof_rotatePnt(value.coord.x,value.coord.y,theta)
        value.coord.x=xy[1]
        value.coord.y=xy[2]
        if value.coord.x<os
            os=value.coord.x
        end
    end
    return os
end

#applies rotational matrix individual coordinates
function lof_rotatePnt(x,y,theta)
    co_od=[x y]
    rotated=co_od*[cos(theta) -1*sin(theta);sin(theta) cos(theta)]
    return rotated
end

#translates the entire region by specified offset
#sets unique IDs for owpps and pccs
function lof_slideAxis(ocn,os)
    for value in ocn.gens
        value.coord.x=value.coord.x-os
        value.id="1"*string(value.num)
    end
    for value in ocn.pccs
        value.coord.x=value.coord.x-os
        value.id="2"*string(value.num)
    end
end

#changes angle to an arc length
function lof_deg2lgth(d,dPl)
    return d*dPl
end

#calculates length of 1 deg of longitude at given lattitude
function lof_lg1deg(lat,lngth)
    return cos(lof_d2r(lat))*lngth
end

#finds minimum distance between any 2 owpp
function lof_mnGap(gens)
    lnths=Array{Float64,1}()
    for gen0 in gens
        for gen1 in gens
            if gen0.id != gen1.id
                push!(lnths,lof_pnt2pnt_dist(gen0.coord,gen1.coord))
            end
        end
    end
    mnGp=findmin(lnths)[1]
    return mnGp
end
################################################################################
############################ Cartesian to GPS transform ########################
################################################################################
#location=ocean.gens
#value=location[1]
#ocn=ocean
#lof_unXformAxis(ocean)
#finds original gps corordinates from untransformed cartesian
function lof_cartesian2gps(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.gps.lat=lof_lgth2deg(value.coord.y,lnthLT)+base.lat
        lnthLG=lof_lg1deg(value.gps.lat,lnthLT)
        value.gps.lng=lof_lgth2deg(value.coord.x,lnthLG)+base.lng
    end
end

#changes arc length to angle
function lof_lgth2deg(d,dPl)
    return d/dPl
end

#performs inverse transforms on cartesian coordinates
function lof_unXformAxis(ocn)
    os=ocn.offset
    lof_unSlideAxis(ocn,os)
    lof_unRotateAxis(ocn)
end

#translates the oss by specified offset
function lof_unSlideAxis(ocn,os)
    #for value in ocn.osss
    #    value.coord.x=value.coord.x+os
    #end
    for value in ocn.gens
        value.coord.x=value.coord.x+os
    end
    for value in ocn.pccs
        value.coord.x=value.coord.x+os
    end
end
#lof_unRotateAxis(ocean)
#inverse rotation of oss
function lof_unRotateAxis(ocn)
    rads=(-1)*ocn.angle
    #offset=lof_rotateGroup(ocn.osss,rads,0.0)
    offset=lof_rotateGroup(ocn.gens,rads,0.0)
    offset=lof_rotateGroup(ocn.pccs,rads,0.0)
end
##############################################################################
############################ laying PCC nodes ################################
##############################################################################
#Places the pccs
function lof_layPccs()
    pccs=lod_pccGps()
    location=Array{node,1}()
    for (index, value) in enumerate(pccs)
        shore=node()
        shore.gps.lng=value[1]
        shore.gps.lat=value[2]
        shore.kv=lod_pccKv()
        shore.num=index
        push!(location,deepcopy(shore))
    end
    return location
end
##############################################################################
############################ laying OWPP nodes ################################
##############################################################################
function lof_layGens(ocn)
    num=length(ocn.pccs)
    gpss,mvas,wnds=lod_gensGps()
    locations=Array{node,1}()
    for i=1:length(gpss)
        concession=node()
        concession.gps.lng=gpss[i][1]#set longitude
        concession.gps.lat=gpss[i][2]#set latittude
        concession.mva=mvas[i]#set concession power
        concession.name=wnds[i]#set wind profile name
        concession.kv=lod_cncsKv()#set collector kv
        concession.num=num+i
        push!(locations,deepcopy(concession))
    end
#sorts the owpps into closest to furthest from PCCs
    locations=lof_srtNear2Far(ocn.pccs,locations)
    return locations
end

#sorts the owpps from closest to furthest
function lof_srtNear2Far(pccs,gens)
    #creates a tuple of distances and gen/pcc numbers
    ds=Array{Tuple,1}()
    for gn in gens
        bsf=lof_gps2gps_dist(gn.gps,pccs[1].gps)
        gnn=gn.num
        pcn=pccs[1].num
        for pc in pccs
            if lof_gps2gps_dist(gn.gps,pc.gps)<bsf
                bsf=lof_gps2gps_dist(gn.gps,pc.gps)
                gnn=deepcopy(gn.num)
                pcn=deepcopy(pc.num)
            end
        end
        push!(ds,(bsf,gnn,pcn))
    end

    #sorts tuple by the length entry
    lnths = [x[1] for x in ds]
    ordrd=Array{Tuple,1}()
    for i=1:length(ds)
        ind=findmin(lnths)[2]
        lnths[ind]=Inf
        push!(ordrd,deepcopy(ds[ind]))
    end

    #sorts gens in same order as lengths
    ogens=Array{node,1}()
    for o in ordrd
        for gn in gens
            if gn.num == o[2]
                push!(ogens,gn)
            end
        end
    end

    #re-numbers each owpp
    for i=1:length(ogens)
        ogens[i].num=deepcopy(i+length(pccs))
    end
    return ogens
end

################################################################################
########################### General purpose ####################################
################################################################################
#Change radians to degrees
function lof_r2d(rad)
    return rad*180/pi
end

#Change degrees to radians
function lof_d2r(deg)
    return deg*pi/180
end

#returns the hypotenuse distance between 2 sets of gps coords
function lof_gps2gps_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.lng-pnt1.lng)^2+(pnt2.lat-pnt1.lat)^2)
    return hyp
end