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

typealias RegionGraph DefaultOrderedDict{Int, Dict{Int, atomic_edge}}
typealias BoundingBoxes Dict{Int, Array{Int, 1}}
typealias SupervoxelSize DefaultOrderedDict{Int,Int}
typealias SemanticInfo DefaultOrderedDict{Int,Array{Float32,1}}

function read_bboxes(fn)
    bbox_file = open(fn)
    bboxes = BoundingBoxes()
    for ln in eachline(bbox_file)
        data = [parse(Int,x) for x in split(ln, " ")]
        bboxes[data[1]] = data[2:end]
    end
    return bboxes
end

function read_size(fn)
    size_file = open(fn)
    d_sizes = SupervoxelSize(0)
    for ln in eachline(size_file)
        data = split(ln, " ")
        seg_id = parse(Int, data[1])
        sz = parse(Int,data[2])
        d_sizes[seg_id] = sz
    end
    close(size_file)
    return d_sizes
end

function find_facesegs(bboxes, chunk_bbox)
    d_facesegs = Set{Int}()
    for k in keys(bboxes)
        bbox = bboxes[k]
        if any(i->(bbox[i] == chunk_bbox[i]), 1:6)
            push!(d_facesegs, k)
        end
    end
    return d_facesegs
end

function read_semantic(fn)
    sem_file = open(fn)
    d_sem = SemanticInfo(()-> Float32[])
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

function read_rg(fn, pd)
    rg_file = open(fn)
    rg_volume = RegionGraph(()->Dict{Int, atomic_edge}())
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
        #p1 = get(pd, u1, u1)
        #p2 = get(pd, u2, u2)
        rg_volume[p1][p2] = a_edge
        rg_volume[p2][p1] = a_edge
    end
    close(rg_file)
    return num_seg,rg_volume
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

