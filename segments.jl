type atomic_edge
    p1::UInt64
    p2::UInt64
    sum::Float64
    num::Float64
    v1::UInt64
    v2::UInt64
    aff::Float64
    area::Float64
end

typealias RegionGraph DefaultDict{UInt64, Dict{UInt, atomic_edge}}
typealias RegionGraphArray Array{Dict{UInt, atomic_edge},1}
typealias BoundingBoxes Array{Array{Int32, 1}}
typealias SupervoxelSizes Array{UInt64}
typealias SemanticInfo Array{Array{Float64,1}}
typealias SegmentDict Dict{UInt64,UInt64}
typealias SegmentArray Array{UInt64,1}
typealias SupervoxelDict Dict{UInt64, Set{UInt64}}
typealias SupervoxelSet Set{UInt64}


type SupervoxelInfo
    totalSupervoxels::UInt64
    regionGraph::RegionGraphArray
    boundingBoxes::BoundingBoxes
    supervoxelSizes::SupervoxelSizes
    semanticInfo::SemanticInfo
    segmentDict::SegmentArray
    boundarySupervoxels::SupervoxelSet
end

type SegmentInfo
    regionGraph::RegionGraph
    supervoxelDict::SupervoxelDict
end

function read_size_line(ln)
    ss = split(ln)
    return parse(UInt64, ss[1]), parse(UInt64, ss[2])
end

function read_bbox_line(ln)
    ss = split(ln)
    return parse(UInt64, ss[1]), [parse(Int32, x) for x in ss[2:7]]
end

function read_sem_line(ln)
    ss = split(ln)
    return parse(UInt64, ss[1]), [parse(Float64, x) for x in ss[2:5]]
end

function read_rg_line(ln)
    types = (UInt64, UInt64, Float64, Float64, UInt64, UInt64, Float64, Float64)
    ss = split(ln)
    return atomic_edge(ntuple(i->parse(types[i],ss[i]),8)...)
end

function read_bboxes(fn,max_seg)
    bbox_file = open(fn)
    bboxes = BoundingBoxes(max_seg)
    lines = readlines(bbox_file)
    close(bbox_file)
    num_line = length(lines)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid,bbox_array = read_bbox_line(lines[i])
        bboxes[segid] = bbox_array
    end
    lines = []
    gc()
    return bboxes
end

function read_size(fn,max_seg)
    size_file = open(fn)
    d_sizes = SupervoxelSizes(max_seg)
    lines = readlines(size_file)
    close(size_file)
    num_line = length(lines)
    segids = zeros(UInt64,max_seg)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid, size = read_size_line(lines[i])
        d_sizes[segid] = size
        segids[i] = segid
    end
    lines = []
    gc()
    return segids, d_sizes
end

function find_facesegs(segids, bboxes, chunk_bbox)
    d_facesegs = SupervoxelSet()
    for k in segids
        if k == 0
            continue
        end
        bbox = bboxes[k]
        if any(i->(bbox[i] == chunk_bbox[i]), 1:6)
            push!(d_facesegs, k)
        end
    end
    return d_facesegs
end

function read_semantic(fn, max_seg)
    sem_file = open(fn)
    d_sem = SemanticInfo(max_seg)
    lines = readlines(sem_file)
    close(sem_file)
    num_line = length(lines)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid, sem_array = read_sem_line(lines[i])
        d_sem[segid] = sem_array
    end
    lines = []
    gc()
    return d_sem
end

function read_rg_array(fn)
    rg_file = open(fn)
    lines = readlines(rg_file)
    close(rg_file)
    data = split(lines[1])
    num_seg = parse(UInt64, data[2])
    num_edge = parse(UInt64, data[3])
    rg_volume = RegionGraphArray(num_seg)
    for i in 1:num_seg
        rg_volume[i] = Dict{UInt64, atomic_edge}()
    end
    println("Number of threads = $(nthreads())")
    edges = Array{atomic_edge}(num_edge)
    @threads for i in 1:num_edge
        if i % 100000 == 0
            println("reading $i lines")
        end
        edges[i] = read_rg_line(lines[i+1])
    end
    for i in 1:num_edge
        if i % 100000 == 0
            println("generating $i lines")
        end
        a_edge = edges[i]
        rg_volume[a_edge.p1][a_edge.p2] = a_edge
        rg_volume[a_edge.p2][a_edge.p1] = a_edge
    end
    edges = []
    lines = []
    gc()
    return num_seg,rg_volume
end

function read_rg(fn)
    rg_file = open(fn)
    rg_volume = RegionGraph(()->Dict{UInt64, atomic_edge}())
    lines = readlines(rg_file)
    close(rg_file)
    data = split(lines[1])
    num_seg = parse(UInt64, data[2])
    num_edge = parse(UInt64, data[3])
    sizehint!(rg_volume, num_seg)
    println("Number of threads = $(nthreads())")
    edges = Array{atomic_edge}(num_edge)
    @threads for i in 1:num_edge
        if i % 100000 == 0
            println("reading $i lines")
        end
        edges[i] = read_rg_line(lines[i+1])
    end
    for i in 1:num_edge
        if i % 100000 == 0
            println("generating $i lines")
        end
        a_edge = edges[i]
        rg_volume[a_edge.p1][a_edge.p2] = a_edge
        rg_volume[a_edge.p2][a_edge.p1] = a_edge
    end
    edges = []
    lines = []
    gc()
    return num_seg,rg_volume
end

function mst_entry(l)
    data = split(l)
    parent = parse(UInt64, data[3])
    child1 = parse(UInt64, data[1])
    child2 = parse(UInt64, data[2])
    if parent == child1
        return parent, child2
    else
        return parent, child1
    end
end


function agglomerate(fn, max_seg)
    pd = SegmentArray(max_seg)
    segs = SupervoxelDict()
    mst = open(fn)
    lines = readlines(mst)
    num_line = length(lines)
    sizehint!(segs, max_seg-num_line)
    close(mst)
    @threads for i in 1:max_seg
        pd[i] = i
    end
    @threads for i in 1:num_line
        entry = mst_entry(lines[i])
        # the first one is child, the second one is parent
        pd[entry[2]] = entry[1]
    end

    # find the root id
    for c in 1:max_seg
        p = pd[c]
        if p == c
            continue
        end
        # list of child node, for path compression
        clst = Set{UInt64}(c)
        # find the root
        while pd[p] != p
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
@time    max_seg, rg_volume = read_rg_array("rg_volume.in")
@time    _, new_rg = read_rg("new_rg.in")
@time    segs, pd = agglomerate("test_mst.in", max_seg)
@time    segids, d_sizes = read_size("sv_volume.in",max_seg)
@time    d_sem = read_semantic("sem_volume.in",max_seg)
    #rg_faces, face_segs, d_facesegs = process_faces(sgm.segmentation)
@time    bboxes = read_bboxes("bbox_volume.in",max_seg)
    #println("$(length(keys(d_sizes)))")
    #println("size of rg: $(length(keys(new_rg)))")
    #d_facesegs = DefaultOrderedDict{Int,Int}(0)
    println(length(segids))
@time    d_facesegs = find_facesegs(segids, bboxes,chunk_bbox)
    println(length(d_facesegs))
    svInfo = SupervoxelInfo(max_seg, rg_volume, bboxes, d_sizes, d_sem, pd, d_facesegs)
    segInfo = SegmentInfo(new_rg, segs)
    return svInfo, segInfo
end
