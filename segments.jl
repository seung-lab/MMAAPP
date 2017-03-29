type atomic_edge
    p1::Int
    p2::Int
    sum::Float64
    num::Float64
    v1::Int
    v2::Int
    aff::Float64
    area::Float64
end

function read_bboxes(fn)
    bbox_file = open(fn)
    bboxes = Dict{Int, Array{Int, 1}}()
    for ln in eachline(bbox_file)
        data = [parse(Int,x) for x in split(ln, " ")]
        bboxes[data[1]] = data[2:end]
    end
    return bboxes
end

function process_size()
    size_file = open("sv_size.in")
    d_sizes = DefaultOrderedDict{Int,Int}(()->0)
    for ln in eachline(size_file)
        data = split(ln, " ")
        seg_id = parse(Int, data[1])
        sz = parse(Int,data[2])
        d_sizes[seg_id] = sz
    end
    close(size_file)
    return d_sizes
end

function process_faces(seg)
    sz = size(seg)
    rg_faces = DefaultOrderedDict{Int,Set{Int}}(()->Set{Int}())
    d_faceareas = DefaultOrderedDict{Int,Int}(0)

    face_segs = Dict("xy-"=>Set{Int}(), "xy+"=>Set{Int}(), "xz-"=>Set{Int}(), "xz+"=>Set{Int}(), "yz-"=>Set{Int}(), "yz+"=>Set{Int}())

    for z in [1, sz[3]]
        for y in 1:sz[2]
            for x in 1:sz[1]
                d_faceareas[seg[x,y,z]] += 1
                if ( (x > 1) && seg[x-1,y,z]!=0 && seg[x,y,z]!=seg[x-1,y,z])
                  p = minmax(seg[x,y,z], seg[x-1,y,z])
                  push!(rg_faces[p[1]],p[2])
                  push!(rg_faces[p[2]],p[1])
                end
                if ( (y > 1) && seg[x,y-1,z]!=0 && seg[x,y,z]!=seg[x,y-1,z])
                  p = minmax(seg[x,y,z], seg[x,y-1,z])
                  push!(rg_faces[p[1]],p[2])
                  push!(rg_faces[p[2]],p[1])
                end
                if z == 1
                    push!(face_segs["xy-"], seg[x,y,z])
                else
                    push!(face_segs["xy+"], seg[x,y,z])
                end
            end
        end
    end

    for z in 1:sz[3]
        for y in [1,sz[2]]
            for x in 1:sz[1]
                d_faceareas[seg[x,y,z]] += 1
                if ( (x > 1) && seg[x-1,y,z]!=0 && seg[x,y,z]!=seg[x-1,y,z])
                  p = minmax(seg[x,y,z], seg[x-1,y,z])
                  push!(rg_faces[p[1]],p[2])
                  push!(rg_faces[p[2]],p[1])
                end
                if ( (z > 1) && seg[x,y,z-1]!=0 && seg[x,y,z]!=seg[x,y,z-1])
                  p = minmax(seg[x,y,z], seg[x,y,z-1])
                  push!(rg_faces[p[1]],p[2])
                  push!(rg_faces[p[2]],p[1])
                end
                if y == 1
                    push!(face_segs["xz-"], seg[x,y,z])
                else
                    push!(face_segs["xz+"], seg[x,y,z])
                end
            end
        end
    end

    for z in 1:sz[3]
        for y in 1:sz[2]
            for x in [1,sz[1]]
                d_faceareas[seg[x,y,z]] += 1
                if ( (y > 1) && seg[x,y-1,z]!=0 && seg[x,y,z]!=seg[x,y-1,z])
                  p = minmax(seg[x,y,z], seg[x,y-1,z])
                  push!(rg_faces[p[1]],p[2])
                  push!(rg_faces[p[2]],p[1])
                end
                if ( (z > 1) && seg[x,y,z-1]!=0 && seg[x,y,z]!=seg[x,y,z-1])
                  p = minmax(seg[x,y,z], seg[x,y,z-1])
                  push!(rg_faces[p[1]],p[2])
                  push!(rg_faces[p[2]],p[1])
                end
                if x == 1
                    push!(face_segs["yz-"], seg[x,y,z])
                else
                    push!(face_segs["yz+"], seg[x,y,z])
                end
            end
        end
    end
    return rg_faces, face_segs, d_faceareas
end

function process_semantic()
    sem_file = open("sem_volume.in")
    d_sem = DefaultOrderedDict{Int,Array{Float32,1}}(()-> Float32[])
    for ln in eachline(sem_file)
        data = split(ln, " ")
        seg_id = parse(Int, data[1])
        axon = parse(Float32,data[2])
        dend = parse(Float32,data[3])
        glial = parse(Float32,data[4])
        psd = parse(Float32,data[5])
        d_sem[seg_id]=[axon,dend,glial,psd]
    end
    return d_sem
end

function process_volume()
    rg_file = open("rg_volume.in")
    rg_volume = DefaultOrderedDict{Int, Dict{Int, atomic_edge}}(()->Dict{Int, atomic_edge}())
    num_seg = 0
    for ln in eachline(rg_file)
        data = split(ln, " ")
        if length(data) == 3
            num_seg = parse(Int, data[2])
            continue
        end
        u1 = parse(Int, data[5])
        u2 = parse(Int, data[6])
        aff = parse(Float64, data[7])
        area = parse(Float64, data[8])
        s = parse(Float64, data[3])
        n = parse(Float64, data[4])
        p1 = parse(Int, data[1])
        p2 = parse(Int, data[2])
        a_edge = atomic_edge(p1,p2,s,n,u1,u2,aff,area)
        rg_volume[p1][p2] = a_edge
        rg_volume[p2][p1] = a_edge
    end
    close(rg_file)
    return num_seg,rg_volume
end


function read_rg(fn, pd)
    rg_file = open(fn)
    #rg = Dict{Tuple{Int,Int},atomic_edge}()
    rg = DefaultOrderedDict{Int, Dict{Int, atomic_edge}}(()->Dict{Int, atomic_edge}())
    for ln in eachline(rg_file)
        data = split(ln, " ")
        u1 = parse(Int, data[5])
        u2 = parse(Int, data[6])
        aff = parse(Float64, data[7])
        area = parse(Float64, data[8])
        s = parse(Float64, data[3])
        n = parse(Float64, data[4])
        p1 = parse(Int, data[1])
        p2 = parse(Int, data[2])
        a_edge = atomic_edge(p1,p2,s,n,u1,u2,aff,area)
        #p1 = get(pd, u1, u1)
        #p2 = get(pd, u2, u2)
        rg[p1][p2] = a_edge
        rg[p2][p1] = a_edge
    end
    return rg
end

function agglomerate(sgm)
    pd = Dict{UInt32,UInt32}()
    segs = Dict{Int, Set{Int}}()
    for idx in 1:length(sgm.segmentPairAffinities)
        if sgm.segmentPairAffinities[idx] > agg_threshold
            # the first one is child, the second one is parent
            pd[sgm.segmentPairs[idx,1]] = sgm.segmentPairs[idx,2]
        end
    end

    # find the root id
    for (c,p) in pd
        # list of child node, for path compression
        clst = Set{Int}(c)
        # find the root
        while haskey(pd, p)
            push!(clst, p)
            p = pd[p]
        end
        for c in clst
            pd[c] = p
        end
        # now p is the root id
        # path compression
        if haskey(segs, p)
            #println("root already exist")
            union!(segs[p], clst)
        else
            segs[p] = clst
            push!(segs[p],p)
        end
    end
    return segs, pd
end

