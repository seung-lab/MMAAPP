using HDF5
using EMIRT
using DataStructures

type atomic_edge
    v1::Int
    v2::Int
    sum::AbstractFloat
    num::Int
end

function check_edge(edge, seg1, seg2)
    neighboor1 = rg_volume[edge.v1]
    if length(collect(intersect(seg2,neighboor1))) > 1
        return false
    end
    neighboor2 = rg_volume[edge.v2]
    if length(collect(intersect(seg1,neighboor2))) > 1
        return false
    end
    if edge.num > 300 && edge.sum/edge.num < 0.15
        return false
    end
    return true
end

function read_rg(fn)
    rg_file = open(fn)
    rg = Dict{Tuple{Int,Int},atomic_edge}()
    for ln in eachline(rg_file)
        data = split(ln, " ")
        u1 = parse(Int, data[5])
        u2 = parse(Int, data[6])
        s = parse(Float32, data[3])
        n = parse(Int, data[4])
        a_edge = atomic_edge(u1,u2,s,n)
        p1 = parse(Int, data[1])
        p2 = parse(Int, data[2])
        rg[(p1,p2)] = a_edge
    end
    return rg
end


function sum_vol(clst, d_sizes)
    total_vol = 0
    for c in clst
        total_vol += d_sizes[c]
    end
    return total_vol
end

function faces_touched(segs, face_segs)
    faces = Set()
    for k in keys(face_segs)
        if !isempty(intersect(segs,face_segs[k]))
            push!(faces, k)
        end
    end
    println("faces: $(faces)")
    return faces
end

function check_connectivity(s, segment, rg_volume)
    cc = 0
    visited = s
    for root in segment
        if root in visited
            continue
        end
        queue = Queue(Int)
        enqueue!(queue, root)

        while length(queue) > 0
            root = dequeue!(queue)
            if root in visited
                continue
            end

            push!(visited,root)
            for neighboor in rg_volume[root]
                if !(neighboor in visited) && neighboor in segment
                    enqueue!(queue, neighboor)
                end
            end
        end
        cc += 1
    end
    return cc
end

function really_check_freeends(ends, segment, rg_volume, d_sizes)
    free_ends = Set{Int}()
    avg_volume = sum_vol(segment,d_sizes)/length(collect(segment))
    for s in ends
        if d_sizes[s] > 5*avg_volume
            continue
        end
        end_segment = Set{Int}()
        for neighboor in rg_volume[s]
            if neighboor in segment
                push!(end_segment, neighboor)
            end
        end
        push!(end_segment, s)
        cc = check_connectivity(end_segment, segment, rg_volume)
        if cc == 1
            push!(free_ends, s)
        end
    end
    return free_ends
end

function check_segment(segment, rg_volume, d_sizes, d_faceareas)
    free_ends = Set{Int}()
    count = 0
    total_vol = sum_vol(segment, d_sizes)
    for s in segment
        cc = check_connectivity(Set{Int}(s), segment, rg_volume)
        if cc > 10 || (cc > 5 && d_sizes[s] > 0.5*total_vol)
            return Set{Int}()
        end
        if cc == 2
            count += 1
        elseif cc == 1
            push!(free_ends, s)
        end
    end
    if count * 2 > length(collect(segment)) && length(collect(segment)) > 5
        ends = setdiff(free_ends, keys(d_faceareas))
        return really_check_freeends(ends, segment, rg_volume, d_sizes)
    elseif count * 3 > length(collect(segment)) && length(collect(segment)) > 10
        ends = setdiff(free_ends, keys(d_faceareas))
        return really_check_freeends(ends, segment, rg_volume, d_sizes)
    else
        return Set{Int}()
    end
end

function connected(subset, rg, facesizes)
    visited = Set{Int}()
    if length(collect(subset)) == 0
        println("no exit face")
        return Int[]
    end

    ccsz = Int[]
    for root in subset
        if root in visited
            continue
        end
        seg_sz = 0
        queue = Queue(Int)
        enqueue!(queue, root)

        while length(queue) > 0
            root = dequeue!(queue)
            if root in visited
                continue
            end
            seg_sz += facesizes[root]

            push!(visited,root)
            for neighboor in rg[root]
                if !(neighboor in visited) && neighboor in subset
                    enqueue!(queue, neighboor)
                end
            end
        end
        push!(ccsz, seg_sz)
    end
    println("connected components: $ccsz")
    return ccsz
end

function process_volume(seg)
    rg_file = open("rg_volume.in")
    rg_volume = DefaultDict(Int,Set{Int}, ()->Set{Int}())
    for ln in eachline(rg_file)
        data = split(ln, " ")
        u1 = parse(Int,data[1])
        u2 = parse(Int,data[2])
        push!(rg_volume[u1],u2)
        push!(rg_volume[u2],u1)
    end
    close(rg_file)
    size_file = open("sv_size.in")
    d_sizes = DefaultDict(Int,Int,()->0)
    for ln in eachline(size_file)
        data = split(ln, " ")
        seg_id = parse(Int, data[1])
        sz = parse(Int,data[2])
        d_sizes[seg_id] = sz
    end
    close(size_file)
    return rg_volume, d_sizes
end

function process_faces(seg)
    sz = size(seg)
    rg_faces = DefaultDict(Int,Set{Int}, ()->Set{Int}())
    d_faceareas = DefaultDict(Int,Int,0)

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

function agglomerate(sgm)
    thd = 0.2
    pd = Dict{UInt32,UInt32}()
    segs = Dict{Int, Set{Int}}()
    for idx in 1:length(sgm.segmentPairAffinities)
        if sgm.segmentPairAffinities[idx] > thd
            # the first one is child, the second one is parent
            pd[sgm.segmentPairs[idx,1]] = sgm.segmentPairs[idx,2]
        end
    end

    # find the root id
    for (c,p) in pd
        # list of child node, for path compression
        clst = Set{Int}([c])
        # find the root
        while haskey(pd, p)
            push!(clst, p)
            p = pd[p]
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
    return segs
end


sgm = readsgm("sgm.h5")
new_rg = read_rg("new_rg.txt")
println("size of rg: $(length(collect(keys(new_rg))))")
@time rg_volume, d_sizes = process_volume(sgm.segmentation)
println("$(length(collect(keys(d_sizes))))")
segs = agglomerate(sgm)

l_segs = []

for k in keys(segs)
    push!(l_segs, [k,length(collect(segs[k])),sum_vol(collect(segs[k]), d_sizes)])
end

#sort!(l_segs, by=x->x[3]/x[2],rev=true)

rg_faces, face_segs, d_faceareas = process_faces(sgm.segmentation)
println("size of the rg: $(length(collect(keys(rg_faces))))")
checks = []
for l in l_segs
    if l[3] < 10000000 && l[2] > 5
        println("segid: $(l[1]), parts: $(l[2]), size: $(l[3])")
        ccsz = connected(intersect(keys(rg_faces), segs[l[1]]), rg_faces, d_faceareas)
        faces = faces_touched(segs[l[1]], face_segs)
        axion_freeends = check_segment(segs[l[1]], rg_volume, d_sizes, d_faceareas)
        if !isempty(axion_freeends)
            push!(checks, [l[1], axion_freeends])
        end
    end
end
println([x[1] for x in checks])
pairs = []
for p in keys(new_rg)
    a_edge = new_rg[p]
    seg1 = 0
    seg2 = 0
    for i in checks
        if a_edge.v1 in i[2]
            seg1 = i[1]
        elseif a_edge.v2 in i[2]
            seg2 = i[1]
        end
    end
    #if !isempty(seg1) && !isempty(seg2) && isempty(intersect(seg1[2],seg2[2]))
    if seg1 != 0 && seg2 != 0 && check_edge(a_edge, segs[seg1], segs[seg2])
        println("$seg1, $seg2")
        println(a_edge)
        push!(pairs, seg1)
        push!(pairs, seg2)
    end
end
println(pairs)
