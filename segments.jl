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
typealias SupervoxelSizes DefaultOrderedDict{Int,Int}
typealias SemanticInfo DefaultOrderedDict{Int,Array{Float32,1}}
typealias SegmentDict Dict{Int,Int}
typealias SupervoxelDict Dict{Int, Set{Int}}
typealias SupervoxelSet Set{Int}


type SupervoxelInfo
    totalSupervoxels::Int64
    regionGraph::RegionGraph
    boundingBoxes::BoundingBoxes
    supervoxelSizes::SupervoxelSizes
    semanticInfo::SemanticInfo
    segmentDict::SegmentDict
    boundarySupervoxels::SupervoxelSet
end

type SegmentInfo
    regionGraph::RegionGraph
    supervoxelDict::SupervoxelDict
end

function read_bboxes(fn)
    bbox_file = open(fn)
    bboxes = BoundingBoxes()
    lines = readlines(bbox_file)
    for ln in lines
        data = [parse(Int,x) for x in split(ln, " ")]
        bboxes[data[1]] = data[2:end]
    end
    close(bbox_file)
    return bboxes
end

function read_size(fn)
    size_file = open(fn)
    d_sizes = SupervoxelSizes(0)
    lines = readlines(size_file)
    for ln in lines
        data = split(ln)
        seg_id = parse(Int, data[1])
        sz = parse(Int,data[2])
        d_sizes[seg_id] = sz
    end
    close(size_file)
    return d_sizes
end

function find_facesegs(bboxes, chunk_bbox)
    d_facesegs = SupervoxelSet()
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
    lines = readlines(sem_file)
    for ln in lines
        data = split(ln)
        seg_id = parse(Int, data[1])
        axon = parse(Float32,data[2])
        dend = parse(Float32,data[3])
        glial = parse(Float32,data[4])
        psd = parse(Float32,data[5])
        d_sem[seg_id]=[axon,dend,glial,psd]
    end
    close(sem_file)
    return d_sem
end

function read_rg_line(ln)
    types = (Int, Int, Float64, Float64, Int, Int, Float64, Float64)
    ss = split(ln)
    return atomic_edge(ntuple(i->parse(types[i],ss[i]),8)...)
end

function read_rg(fn, pd)
    rg_file = open(fn)
    rg_volume = RegionGraph(()->Dict{Int, atomic_edge}())
    num_seg = 0
    lines = readlines(rg_file)
    data = split(lines[1])
    num_seg = parse(Int, data[2])
    for ln in lines[2:end]
        a_edge = read_rg_line(ln)
        #u1 = parse(Int, data[5])
        #u2 = parse(Int, data[6])
        #aff = parse(Float64, data[7])
        #area = parse(Float64, data[8])
        #s = parse(Float64, data[3])
        #n = parse(Float64, data[4])
        #p1 = parse(Int, data[1])
        #p2 = parse(Int, data[2])
        #a_edge = atomic_edge(p1,p2,s,n,u1,u2,aff,area)
        #p1 = get(pd, u1, u1)
        #p2 = get(pd, u2, u2)
        rg_volume[a_edge.p1][a_edge.p2] = a_edge
        rg_volume[a_edge.p2][a_edge.p1] = a_edge
    end
    close(rg_file)
    return num_seg,rg_volume
end

function mst_entry(l)
    data = split(l)
    parent = parse(Int, data[3])
    child1 = parse(Int, data[1])
    child2 = parse(Int, data[2])
    weight = parse(Float64, data[4])/parse(Float64, data[5])
    if parent == child1
        return (parent, child2, weight)
    else
        return (parent, child1, weight)
    end
end


function agglomerate(fn)
    pd = SegmentDict()
    segs = SupervoxelDict()
    mst = open(fn)
    lines = readlines(mst)
    for l in lines
        entry = mst_entry(l)
            # the first one is child, the second one is parent
            pd[entry[2]] = entry[1]
    end
    close(mst)

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

function load_segments(chunk_bbox)
    segs, pd = agglomerate("test_mst.in")
    num_seg, rg_volume = read_rg("rg_volume.in", Dict{UInt32,UInt32}())
    d_sizes = read_size("sv_size.in")
    d_sem = read_semantic("sem_volume.in")
    #rg_faces, face_segs, d_facesegs = process_faces(sgm.segmentation)
    _, new_rg = read_rg("new_rg.in", pd)
    bboxes = read_bboxes("bbox_volume.in")
    println("$(length(keys(d_sizes)))")
    println("size of rg: $(length(keys(new_rg)))")
    #d_facesegs = DefaultOrderedDict{Int,Int}(0)
    d_facesegs = find_facesegs(bboxes,chunk_bbox)
    svInfo = SupervoxelInfo(num_seg, rg_volume, bboxes, d_sizes, d_sem, pd, d_facesegs)
    segInfo = SegmentInfo(new_rg, segs)
    return svInfo, segInfo
end
