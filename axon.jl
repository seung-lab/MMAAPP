using HDF5
using EMIRT
using DataStructures
using Combinatorics

type atomic_edge
    p1::Int
    p2::Int
    sum::Float32
    num::Int
    v1::Int
    v2::Int
    aff::Float32
    area::Int
end

function print_edge(edge, mean_aff)
    area = edge.num
    sum_aff = edge.sum
    if mean_aff > 0
        sum_aff = area*mean_aff
    end
    println("$(edge.p1) $(edge.p2) $(sum_aff) $(area) $(edge.v1) $(edge.v2) $(edge.aff) $(edge.area)")
end

function check_edge(edge, seg1, seg2)
    neighboor1 = Set{Int}()
    neighboor2 = Set{Int}()
    if edge.v1 in seg1
        neighboor1 = rg_volume[edge.v1]
        neighboor2 = rg_volume[edge.v2]
    else
        neighboor1 = rg_volume[edge.v2]
        neighboor2 = rg_volume[edge.v1]
    end
    #if length(intersect(seg2,neighboor1)) > 1
    #    return false
    #end
    #if length(intersect(seg1,neighboor2)) > 1
    #    return false
    #end
    #if edge.num > 300 && edge.sum/edge.num < 0.15
    #    return false
    #end
    len1 = length(intersect(seg2,neighboor1))
    len2 = length(intersect(seg1,neighboor2))
    if len1 == 1 && len2 == 1
        if edge.num > 300 && edge.sum/edge.num < 0.15
            return false
        end
    elseif len1 > 1 && len2 > 1
        return false
    else
        if edge.sum/edge.num < 0.10
            return false
        end
    end
    return true
end

function read_rg(fn, pd)
    rg_file = open(fn)
    #rg = Dict{Tuple{Int,Int},atomic_edge}()
    rg = DefaultDict(Int, Dict{Int, atomic_edge}, ()->Dict{Int, atomic_edge}())
    for ln in eachline(rg_file)
        data = split(ln, " ")
        u1 = parse(Int, data[5])
        u2 = parse(Int, data[6])
        aff = parse(Float32, data[7])
        area = parse(Int, data[8])
        s = parse(Float32, data[3])
        n = parse(Int, data[4])
        p1 = parse(Int, data[1])
        p2 = parse(Int, data[2])
        a_edge = atomic_edge(p1,p2,s,n,u1,u2,aff,area)
        p1 = get(pd, u1, u1)
        p2 = get(pd, u2, u2)
        rg[p1][p2] = a_edge
        rg[p2][p1] = a_edge
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

function really_check_freeends(ends, segment, rg_volume, d_sizes, d_faceareas)
    free_ends = Set{Int}()
    boundary_segs = keys(d_faceareas)
    for s in ends
        end_segment = Set{Int}()
        count = 0
        if s in boundary_segs
            continue
        end
        for neighboor in rg_volume[s]
            if neighboor in segment
                count += 1
                push!(end_segment, neighboor)
                if neighboor in boundary_segs
                    count = 100
                    break
                end
            end
        end
        if count > 2
            continue
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
    #ends = setdiff(free_ends, keys(d_faceareas))
    ends = really_check_freeends(free_ends, segment, rg_volume, d_sizes, d_faceareas)
    if length(ends) > 3 && (count * 3 < length(segment))
        return Set{Int}()
    end
    return ends
end

function connected(subset, rg, facesizes)
    visited = Set{Int}()
    if length(subset) == 0
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

function match_axons(axons, segs, new_rg, free_ends)
    visited = Set{atomic_edge}()
    pairs = Int[]
    edges = Dict{atomic_edge, Float32}()
    for a in keys(axons)
        matches = intersect(keys(new_rg[a]),keys(axons))
        #println("test: $a")
        for b in matches
            #println("test: $a, $b")
            a_edge = new_rg[a][b]
            if a_edge in visited
                continue
            end
            push!(visited, a_edge)
            if !(a_edge.v1 in free_ends) || !(a_edge.v2 in free_ends)
                continue
            end
            if check_edge(a_edge, segs[a], segs[b])
                p = minmax(a,b)
                #println("$(p[1]), $(p[2])")
                println(a_edge)
                edges[a_edge] = 0.199999
                push!(pairs, a)
                push!(pairs, b)
            end
        end
    end
    println(pairs)
    return edges
end

function match_long_axons(small_pieces, long_axons, segs, new_rg, free_ends)
    for s in small_pieces
        neighboors = keys(new_rg[s])
        candidates = intersect(neighboors, keys(long_axons))
        axons = []
        for c in candidates
            a_edge = new_rg[s][c]
            if !(a_edge.v1 in free_ends) && !(a_edge.v2 in free_ends)
                continue
            end
            push!(axons, c)
        end
        if length(axons) > 1
            for p in combinations(axons,2)
                if haskey(new_rg[p[1]],p[2])
                    continue
                end
                push!(p, s)
                println("axons: $s $p")
            end
        end
    end
end



sgm = readsgm("sgm.h5")
@time rg_volume, d_sizes = process_volume(sgm.segmentation)
println("$(length(keys(d_sizes)))")
segs, pd = agglomerate(sgm)
new_rg = read_rg("new_rg_2.txt", pd)
println("size of rg: $(length(keys(new_rg)))")

l_segs = []

for k in keys(segs)
    push!(l_segs, [k,length(segs[k]),sum_vol(segs[k], d_sizes)])
end

#sort!(l_segs, by=x->x[3]/x[2],rev=true)

@time rg_faces, face_segs, d_faceareas = process_faces(sgm.segmentation)
println("size of rg: $(length(keys(rg_faces)))")
axons = Dict{Int, Set{Int}}()
long_axons = Dict{Int, Set{Int}}()
free_ends = Set{Int}()
small_pieces = Set{Int}()
for l in l_segs
    if l[3] < 10000000 && l[2] > 5
        #ccsz = connected(intersect(keys(rg_faces), segs[l[1]]), rg_faces, d_faceareas)
        #faces = faces_touched(segs[l[1]], face_segs)
        axon_freeends = check_segment(segs[l[1]], rg_volume, d_sizes, d_faceareas)
        if !isempty(axon_freeends) && length(axon_freeends) < 10
            #println("segid: $(l[1]), parts: $(l[2]), size: $(l[3]), free_ends: $(length(axon_freeends)) ($(axon_freeends))")
            #push!(checks, [l[1], axon_freeends])
            axons[l[1]] = axon_freeends
            union!(free_ends, axon_freeends)
            if l[2] > 20
                long_axons[l[1]] = axon_freeends
            end
        end
    elseif l[2] < 5
        if isempty(intersect(segs[l[1]], keys(d_faceareas)))
            push!(small_pieces, l[1])
        end
    end
end
println("$(length(keys(axons))) axon candidates")
println("$(length(small_pieces)) small pieces")
println("$(length(keys(long_axons))) long axon candidates")
#match_long_axons(small_pieces, long_axons, segs, new_rg, free_ends)
matches = match_axons(axons, segs, new_rg, free_ends)
visited = Set{atomic_edge}()
for edge in keys(matches)
    print_edge(edge, matches[edge])
    push!(visited, edge)
end
for a in keys(new_rg)
    for b in keys(new_rg[a])
        edge = new_rg[a][b]
        if edge in visited
            continue
        else
            print_edge(edge, -1)
            push!(visited, edge)
        end
    end
end
#println([x[1] for x in checks])
