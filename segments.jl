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

function read_size_line(ln)
    return [parse(Int,x) for x in split(ln)]
end

function read_bbox_size_line(ln)
    return [parse(Int,x) for x in split(ln)]
end

function read_sem_line(ln)
    ss = split(ln)
    return parse(Int, ss[1]), [parse(Float64, x) for x in ss[2:5]]
end

function read_rg_line(ln)
    types = (Int, Int, Float64, Float64, Int, Int, Float64, Float64)
    ss = split(ln)
    return atomic_edge(ntuple(i->parse(types[i],ss[i]),8)...)
end

function read_bboxes(fn)
    bbox_file = open(fn)
    bboxes = BoundingBoxes()
    lines = readlines(bbox_file)
    close(bbox_file)
    num_line = length(lines)
    bbox_array = Array{Array{Int,1}}(num_line)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        bbox_array[i] = read_bbox_size_line(lines[i])
    end
    for i in 1:num_line
        if i % 100000 == 0
            println("generating $i lines")
        end
        bboxes[bbox_array[i][1]] = bbox_array[i][2:7]
    end
    return bboxes
end

function read_size(fn)
    size_file = open(fn)
    d_sizes = SupervoxelSizes(0)
    lines = readlines(size_file)
    close(size_file)
    num_line = length(lines)
    size_array = Array{Array{Int,1}}(num_line)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        size_array[i] = read_bbox_size_line(lines[i])
    end
    for i in 1:num_line
        if i % 100000 == 0
            println("generating $i lines")
        end
        d_sizes[size_array[i][1]] = size_array[i][2]
    end
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
    d_sem = SemanticInfo(()-> Float64[])
    lines = readlines(sem_file)
    close(sem_file)
    num_line = length(lines)
    segid = Array{Int}(num_line)
    sem_array = Array{Array{Float64,1}}(num_line)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid[i], sem_array[i] = read_sem_line(lines[i])
    end
    for i in 1:num_line
        if i % 100000 == 0
            println("generating $i lines")
        end
        d_sem[segid[i]]=sem_array[i]
    end
    return d_sem
end

function read_rg(fn)
    rg_file = open(fn)
    rg_volume = RegionGraph(()->Dict{Int, atomic_edge}())
    lines = readlines(rg_file)
    close(rg_file)
    data = split(lines[1])
    num_seg = parse(Int, data[2])
    num_edge = parse(Int, data[3])
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
    close(mst)
    for l in lines
        entry = mst_entry(l)
        # the first one is child, the second one is parent
        pd[entry[2]] = entry[1]
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

function load_segments(chunk_bbox)
@time    num_seg, rg_volume = read_rg("rg_volume.in")
@time    _, new_rg = read_rg("new_rg.in")
@time    segs, pd = agglomerate("test_mst.in")
@time    d_sizes = read_size("sv_volume.in")
@time    d_sem = read_semantic("sem_volume.in")
    #rg_faces, face_segs, d_facesegs = process_faces(sgm.segmentation)
@time    bboxes = read_bboxes("bbox_volume.in")
    #println("$(length(keys(d_sizes)))")
    #println("size of rg: $(length(keys(new_rg)))")
    #d_facesegs = DefaultOrderedDict{Int,Int}(0)
@time    d_facesegs = find_facesegs(bboxes,chunk_bbox)
    svInfo = SupervoxelInfo(num_seg, rg_volume, bboxes, d_sizes, d_sem, pd, d_facesegs)
    segInfo = SegmentInfo(new_rg, segs)
    return svInfo, segInfo
end
