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
    function SupervoxelInfo()
        totalSupervoxels = 0
        regionGraph = RegionGraphArray()
        boundingBoxes = BoundingBoxes()
        supervoxelSizes = SupervoxelSizes()
        semanticInfo = SemanticInfo()
        segmentDict = SegmentArray()
        boundarySupervoxels = SupervoxelSet()
        return new(totalSupervoxels, regionGraph, boundingBoxes, supervoxelSizes, semanticInfo, segmentDict, boundarySupervoxels)
    end
end

type SegmentInfo
    regionGraph::RegionGraph
    supervoxelDict::SupervoxelDict
    function SegmentInfo()
        regionGraph = RegionGraph(()->Dict{UInt64, atomic_edge}())
        supervoxelDict = SupervoxelDict()
        return new(regionGraph, supervoxelDict)
    end
end

function read_size_line(ln)
    ss = split(ln)
    return parse(UInt64, ss[1]), parse(UInt64, ss[2])
end

function read_bbox_line(ln)
    ss = split(ln)
    return parse(UInt64, ss[1]), [parse(Int32, ss[i]) for i in 2:7]
end

function read_sem_line(ln)
    ss = split(ln)
    return parse(UInt64, ss[1]), [parse(Float64, ss[i]) for i in 2:5]
end

function read_rg_line(ln)
    types = (UInt64, UInt64, Float64, Float64, UInt64, UInt64, Float64, Float64)
    ss = split(ln)
    return atomic_edge(ntuple(i->parse(types[i],ss[i]),8)...)
end

function read_bboxes!(fn,max_seg,svInfo)
    bbox_file = open(fn)
    svInfo.boundingBoxes = BoundingBoxes(max_seg)
    lines = readlines(bbox_file)
    close(bbox_file)
    num_line = length(lines)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid,bbox_array = read_bbox_line(lines[i])
        svInfo.boundingBoxes[segid] = bbox_array
    end
    lines = []
    gc()
end

function read_size!(fn,max_seg,svInfo)
    size_file = open(fn)
    svInfo.supervoxelSizes = SupervoxelSizes(max_seg)
    lines = readlines(size_file)
    close(size_file)
    num_line = length(lines)
    segids = zeros(UInt64,max_seg)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid, size = read_size_line(lines[i])
        svInfo.supervoxelSizes[segid] = size
        segids[i] = segid
    end
    lines = []
    gc()
    return segids
end

function find_facesegs!(segids, chunk_bbox, svInfo)
    svInfo.boundarySupervoxels = SupervoxelSet()
    for k in segids
        if k == 0
            continue
        end
        bbox = svInfo.boundingBoxes[k]
        if any(i->(bbox[i] == chunk_bbox[i]), 1:6)
            push!(svInfo.boundarySupervoxels, k)
        end
    end
end

function read_semantic!(fn, max_seg, svInfo)
    sem_file = open(fn)
    svInfo.semanticInfo = SemanticInfo(max_seg)
    lines = readlines(sem_file)
    close(sem_file)
    num_line = length(lines)
    for i in 1:num_line
        if i % 100000 == 0
            println("reading $i lines")
        end
        segid, sem_array = read_sem_line(lines[i])
        svInfo.semanticInfo[segid] = sem_array
    end
    lines = []
    gc()
end

function read_rg_array!(fn, svInfo)
    rg_file = open(fn)
    lines = readlines(rg_file)
    close(rg_file)
    data = split(lines[1])
    num_seg = parse(UInt64, data[2])
    svInfo.totalSupervoxels = num_seg
    num_edge = parse(UInt64, data[3])
    svInfo.regionGraph = RegionGraphArray(num_seg)
    for i in 1:num_seg
        svInfo.regionGraph[i] = Dict{UInt64, atomic_edge}()
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
        svInfo.regionGraph[a_edge.p1][a_edge.p2] = a_edge
        svInfo.regionGraph[a_edge.p2][a_edge.p1] = a_edge
    end
    edges = []
    lines = []
    gc()
end

function read_rg!(fn, segInfo)
    rg_file = open(fn)
    lines = readlines(rg_file)
    close(rg_file)
    data = split(lines[1])
    num_seg = parse(UInt64, data[2])
    num_edge = parse(UInt64, data[3])
    sizehint!(segInfo.regionGraph, num_seg)
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
        segInfo.regionGraph[a_edge.p1][a_edge.p2] = a_edge
        segInfo.regionGraph[a_edge.p2][a_edge.p1] = a_edge
    end
    edges = []
    lines = []
    gc()
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


function agglomerate!(fn, max_seg, svInfo, segInfo)
    svInfo.segmentDict = SegmentArray(max_seg)
    segInfo.supervoxelDict = SupervoxelDict()
    mst = open(fn)
    lines = readlines(mst)
    num_line = length(lines)
    sizehint!(segInfo.supervoxelDict, max_seg-num_line)
    close(mst)
    @threads for i in 1:max_seg
        svInfo.segmentDict[i] = i
    end
    @threads for i in 1:num_line
        entry = mst_entry(lines[i])
        # the first one is child, the second one is parent
        svInfo.segmentDict[entry[2]] = entry[1]
    end

    # find the root id
    for c in 1:max_seg
        p = svInfo.segmentDict[c]
        if p == c
            continue
        end
        # list of child node, for path compression
        clst = Set{UInt64}(c)
        # find the root
        while svInfo.segmentDict[p] != p
            push!(clst, p)
            p = svInfo.segmentDict[p]
        end
        for c in clst
            svInfo.segmentDict[c] = p
        end
        # now p is the root id
        # path compression
        if haskey(segInfo.supervoxelDict, p)
            #println("root already exist")
            union!(segInfo.supervoxelDict[p], clst)
        else
            segInfo.supervoxelDict[p] = clst
            push!(segInfo.supervoxelDict[p],p)
        end
    end
end

function load_segments!(chunk_bbox, svInfo, segInfo)
@time    read_rg_array!("rg_volume.in", svInfo)
@time    read_rg!("new_rg.in", segInfo)
    max_seg = svInfo.totalSupervoxels
@time    agglomerate!("test_mst.in", max_seg, svInfo, segInfo)
@time    segids = read_size!("sv_volume.in", max_seg, svInfo)
@time    read_semantic!("sem_volume.in", max_seg, svInfo)
    #rg_faces, face_segs, d_facesegs = process_faces(sgm.segmentation)
@time    read_bboxes!("bbox_volume.in", max_seg, svInfo)
    #println("$(length(keys(d_sizes)))")
    #println("size of rg: $(length(keys(new_rg)))")
    #d_facesegs = DefaultOrderedDict{Int,Int}(0)
@time    find_facesegs!(segids, chunk_bbox, svInfo)
end
