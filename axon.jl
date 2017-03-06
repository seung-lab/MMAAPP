using HDF5
using EMIRT
using DataStructures
using Combinatorics

type atomic_edge
    p1::Int
    p2::Int
    sum::Float64
    num::Int
    v1::Int
    v2::Int
    aff::Float64
    area::Int
end

agg_threshold = parse(Float64, ARGS[1])
th_tier1 = agg_threshold - 0.000005
th_tier2 = agg_threshold - 0.000015
th_tier3 = agg_threshold - 0.000025
reliable_th = parse(Float64, ARGS[2])

function print_edge(rg_out, edge, mean_aff)
    area = edge.num
    sum_aff = edge.sum
    if mean_aff > 0
        sum_aff = area*mean_aff
    else
        sum_aff *= 0.9995
    end
    p = minmax(edge.p1, edge.p2)
    push!(rg_out,"$(p[1]) $(p[2]) $(sum_aff) $(area) $(edge.v1) $(edge.v2) $(edge.aff) $(edge.area)\n")
end

function update_new_rg(num_seg, new_rg, matches)
    f = open("axon.in","w")
    visited = Set{atomic_edge}()
    rg_out = AbstractString[]
    for edge in keys(matches)
        print_edge(rg_out, edge, matches[edge])
        push!(visited, edge)
    end
    for a in keys(new_rg)
        for b in keys(new_rg[a])
            edge = new_rg[a][b]
            if edge in visited
                continue
            else
                print_edge(rg_out, edge, -1)
                push!(visited, edge)
            end
        end
    end
    write(f,"$(num_seg-1) $(num_seg) $(length(rg_out))\n")
    for str in rg_out
        write(f, str)
    end
    close(f)
end

function check_edge(edge, seg1, seg2)
    neighboor1 = Set{Int}()
    neighboor2 = Set{Int}()
    if edge.v1 in seg1
        neighboor1 = keys(rg_volume[edge.v1])
        neighboor2 = keys(rg_volume[edge.v2])
    else
        neighboor1 = keys(rg_volume[edge.v2])
        neighboor2 = keys(rg_volume[edge.v1])
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
        if edge.num > 300 && edge.sum/edge.num < reliable_th
            return false
        end
    elseif len1 > 1 && len2 > 1
        return false
    else
        if edge.sum/edge.num < reliable_th
            return false
        end
    end
    return true
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
        area = convert(Int, parse(Float64, data[8]))
        s = parse(Float64, data[3])
        n = convert(Int, parse(Float64, data[4]))
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

function sum_sem(clst, d_sem)
    total_sem = Float32[0,0,0,0]
    for c in clst
        total_sem += d_sem[c]
    end
    return total_sem
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
            for neighboor in keys(rg_volume[root])
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
    avg_vol = sum_vol(segment, d_sizes)/length(segment)
    for s in ends
        end_segment = Set{Int}()
        near_boundary = false
        if s in boundary_segs
            continue
        end
        if d_sizes[s] > 5*avg_vol
            continue
        end
        neighboor2 = Set{Int}()
        for neighboor in keys(rg_volume[s])
            if neighboor in segment
                union!(neighboor2,intersect(keys(rg_volume[neighboor]),segment))
                push!(end_segment, neighboor)
                if neighboor in boundary_segs
                    near_boundary = true
                    break
                end
            end
        end
        if near_boundary || length(neighboor2) > 7
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

function check_dend(segment, rg_volume, d_sizes, l)
    total_vol = l[3]
    sem = l[4]
    if total_vol < 1000000 || sem[2] < sem[1] || sem[2] < 0.5*sem[3]
        return false
    end
    for s in segment
        cc = check_connectivity(Set{Int}(s), segment, rg_volume)
        if (cc > 2 && d_sizes[s] > 0.5*total_vol)
            return true
        end
    end
    return false
end

function check_segment(segment, rg_volume, d_sizes, d_faceareas)
    free_ends = Set{Int}()
    count = 0
    total_vol = sum_vol(segment, d_sizes)
    for s in segment
        cc = check_connectivity(Set{Int}(s), segment, rg_volume)
        if cc > 10 || (cc > 5 && d_sizes[s] > 0.5*total_vol)
            return "dend", Set{Int}()
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
        return "glial", Set{Int}()
    end
    if count * 2 > length(segment)
        return "axon", ends
    else
        return "not sure", ends
    end
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
        area = parse(Int, data[8])
        s = parse(Float64, data[3])
        n = parse(Int, data[4])
        p1 = parse(Int, data[1])
        p2 = parse(Int, data[2])
        a_edge = atomic_edge(p1,p2,s,n,u1,u2,aff,area)
        rg_volume[p1][p2] = a_edge
        rg_volume[p2][p1] = a_edge
    end
    close(rg_file)
    return num_seg,rg_volume
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

function merge_edges(merge_graph, threshold, new_rg)
    edges = OrderedDict{atomic_edge, Float64}()
    visited = Set{Int}()
    for p in keys(merge_graph)
        axons_group = Int[]
        if p in visited
            continue
        end
        queue = Queue(Int)
        enqueue!(queue, p)
        push!(axons_group, p)

        while length(queue) > 0
            root = dequeue!(queue)
            if root in visited
                continue
            end

            push!(visited,root)
            for neighboor in merge_graph[root]
                if !(neighboor in visited)
                    enqueue!(queue, neighboor)
                    push!(axons_group, neighboor)
                end
            end
        end
        for q in combinations(axons_group,2)
            if haskey(new_rg[q[1]], q[2])
                edges[new_rg[q[1]][q[2]]] = threshold
            end
        end
    end
    return edges
end


function match_axons(axons, segs, new_rg, free_ends, considered, is_strict, merge_graph)
    visited = Set{atomic_edge}()
    pairs = Int[]
    processed = Set{Int}()
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
            if !(a_edge.v1 in keys(free_ends)) || a_edge.v1 in considered
                continue
            end
            if !(a_edge.v2 in keys(free_ends)) || a_edge.v2 in considered
                continue
            end
            if is_strict && a_edge.sum/a_edge.num < reliable_th
                continue
            end
            if check_edge(a_edge, segs[a], segs[b])
                p = minmax(a,b)
                println("$(p[1]), $(p[2])")
                println(a_edge)
                push!(processed, a_edge.v1)
                push!(processed, a_edge.v2)
                push!(pairs, a)
                push!(pairs, b)
                push!(merge_graph[a],b)
                push!(merge_graph[b],a)
            end
        end
    end

    println(pairs)
    return processed
end

function match_branches(really_long_axons, long_axons, segs, new_rg, free_ends, considered, merge_graph)
    pairs = []
    processed = Set{Int}()
    for l in really_long_axons
        neighboors = keys(new_rg[l])
        #candidates = intersect(neighboors, really_long_axons)
        candidates = intersect(neighboors, keys(long_axons))
        axons = []
        for c in candidates
            a_edge = new_rg[l][c]
            if (a_edge.v1 in keys(free_ends)) && (a_edge.v2 in keys(free_ends))
                continue
            end
            if (a_edge.v1 in considered) || (a_edge.v2 in considered)
                continue
            end
            if (a_edge.v1 in segs[c]) && !(a_edge.v1 in long_axons[c])
                continue
            end
            if (a_edge.v2 in segs[c]) && !(a_edge.v2 in long_axons[c])
                continue
            end
            if a_edge.sum/a_edge.num > reliable_th
                println("$l, $c")
                push!(pairs,l)
                push!(pairs,c)
                push!(merge_graph[l],c)
                push!(merge_graph[c],l)
                push!(processed, c)
            end
        end
    end
    println(pairs)
    return processed
end

function match_long_axons2(long_axons, new_rg, rg_volume, segs, d_sizes, d_faceareas, considered, merge_graph)
    pairs = []
    processed = Set{Int}()
    for s in keys(d_sizes)
        if s in keys(d_faceareas)
            continue
        end
        if d_sizes[s] > 10000
            continue
        end
        ends = setdiff(intersect(keys(rg_volume[s]), keys(free_ends)),considered)
        if length(ends) > 1
            for p in combinations(ends,2)
                axon1 = free_ends[p[1]]
                axon2 = free_ends[p[2]]
                if axon1 in keys(new_rg[axon2]) || !(axon1 in keys(long_axons)) || !(axon2 in keys(long_axons)) || axon1 == axon2
                    continue
                end
                freeend_dominate = true
                for q in p
                    freeend_area = rg_volume[s][q].area
                    if length(intersect(segs[free_ends[q]], keys(rg_volume[s]))) > 1
                        freeend_dominate = false
                        break
                    end
                    if !freeend_dominate
                        break
                    end
                    #for r in intersect(segs[free_ends[q]], keys(rg_volume[s]))
                    #    q_area = rg_volume[s][r].area
                    #    if q_area > freeend_area
                    #        freeend_dominate = false
                    #        break
                    #    end
                    #end
                end
                if !freeend_dominate
                    continue
                end
                #println("$axon1, $axon2")
                push!(processed, p[1])
                push!(processed, p[2])
                push!(pairs, axon1)
                push!(pairs, axon2)
                push!(merge_graph[axon1], axon2)
                push!(merge_graph[axon2], axon1)
                #update the mean affinity later
                a_edge = atomic_edge(axon1, axon2, 0.2, 1, p[1], p[2], 0.2, 1)
                new_rg[axon1][axon2] = a_edge
                new_rg[axon2][axon1] = a_edge
            end
        end
    end
    println(pairs)
    return processed
end

function match_long_axons(small_pieces, long_axons, new_rg, considered, is_strict, merge_graph)
    processed = Set{Int}()
    pairs = Set{Int}()
    for s in small_pieces
        neighboors = keys(new_rg[s])
        candidates = intersect(neighboors, keys(long_axons))
        axons = []
        for c in candidates
            a_edge = new_rg[s][c]
            free_end = a_edge.v1
            if !(free_end in long_axons[c])
                free_end = a_edge.v2
            end
            if !(free_end in long_axons[c])
                continue
            end
            if free_end in considered
                continue
            end
            push!(axons, [c, free_end, a_edge])
        end
        if length(axons) > 1
            for p in combinations(axons,2)
                if haskey(new_rg[p[1][1]],p[2][1])
                    continue
                end
                if is_strict && p[1][3].sum/p[1][3].num > reliable_th && p[2][3].sum/p[2][3].num > reliable_th
                    push!(pairs, s)
                    push!(pairs, p[1][1])
                    push!(pairs, p[2][1])
                    push!(merge_graph[s],p[1][1])
                    push!(merge_graph[s],p[2][1])
                    push!(merge_graph[p[1][1]], s)
                    push!(merge_graph[p[2][1]], s)
                    push!(processed, p[1][2])
                    push!(processed, p[2][2])
                elseif !is_strict && (p[1][3].sum/p[1][3].num > reliable_th || p[2][3].sum/p[2][3].num > reliable_th)
                    push!(pairs, s)
                    push!(pairs, p[1][1])
                    push!(pairs, p[2][1])
                    push!(merge_graph[s],p[1][1])
                    push!(merge_graph[s],p[2][1])
                    push!(merge_graph[p[1][1]], s)
                    push!(merge_graph[p[2][1]], s)
                    push!(processed, p[1][2])
                    push!(processed, p[2][2])
                end
                println("axons: $s $p")
            end
        end
    end
    println(pairs)
    return processed
end

function process_edge(set_a, set_b, rg_volume)
    for a in set_a
        for b in intersect(keys(rg_volume[a]),set_b)
            edge = rg_volume[a][b]
            if edge.num > 1000
                continue
            end
            if edge.aff/edge.area > 0.5 || (edge.aff/edge.area > agg_threshold && 50 < edge.num < 300)
                println("process: $(length(set_a)), $(length(set_b))")
                println("atomic edge: $edge")
                return true
            end
        end
    end
    return false
end

function process_rg(new_rg, l_segs, segs, rg_volume, d_size, considered, merge_graph)
    visited = Set{atomic_edge}()
    dend_candidates = Set{Int}()
    spine_candidates = Set{Int}()
    queue = Queue(Int)
    for l in l_segs
        if check_dend(segs[l[1]], rg_volume, d_size, l)
            enqueue!(queue, l[1])
        end
    end
    while length(queue) > 0
        a = dequeue!(queue)
        for b in keys(new_rg[a])
            if haskey(segs,b) && (length(segs[b]) > 20 || sum_vol(segs[b], d_size) > 1000000)
                continue
            #elseif !haskey(segs,b) && d_size[b] < 5000
            #    continue
            end
            if a in considered || b in considered
                continue
            end
            edge = new_rg[a][b]
            if edge in visited
                continue
            elseif edge.sum == edge.aff
                continue
            end
            push!(visited, edge)
            set_a = Set([a])
            set_b = Set([b])
            if haskey(segs,a)
                set_a = segs[a]
            end
            if haskey(segs,b)
                set_b = segs[b]
            end
            if process_edge(set_a, set_b, rg_volume) && !haskey(merge_graph,b)
                push!(dend_candidates,a)
                push!(spine_candidates,b)
                push!(merge_graph[a],b)
                push!(merge_graph[b],a)
                enqueue!(queue,b)
            end
        end
    end
    println("dend candidates:")
    println("$dend_candidates, $(length(spine_candidates))")
    return spine_candidates
end



sgm = readsgm("sgm.h5")
@time num_seg, rg_volume = process_volume()
@time d_sizes = process_size()
@time d_sem = process_semantic()
println("$(length(keys(d_sizes)))")
segs, pd = agglomerate(sgm)
new_rg = read_rg("new_rg.in", pd)
println("size of rg: $(length(keys(new_rg)))")

l_segs = []

for k in keys(segs)
    push!(l_segs, [k,length(segs[k]),sum_vol(segs[k], d_sizes),sum_sem(segs[k], d_sem)])
end

#sort!(l_segs, by=x->x[3]/x[2],rev=true)

@time rg_faces, face_segs, d_faceareas = process_faces(sgm.segmentation)
println("size of rg: $(length(keys(rg_faces)))")
axons = Dict{Int, Set{Int}}()
free_ends = Dict{Int, Int}()
small_pieces = Set{Int}()
long_axons = Dict{Int, Set{Int}}()
really_long_axons = Set{Int}()
#glial_segs=Set{Int}()
for l in l_segs
    if l[3] < 10000000 && l[2] > 5
        #ccsz = connected(intersect(keys(rg_faces), segs[l[1]]), rg_faces, d_faceareas)
        #faces = faces_touched(segs[l[1]], face_segs)
        if l[4][3] > l[4][1] && l[4][3] > l[4][2]
            #push!(glial_segs, l[1])
            continue
        end
        seg_type, axon_freeends = check_segment(segs[l[1]], rg_volume, d_sizes, d_faceareas)
        #println("segid: $(l[1]), parts: $(l[2]), size: $(l[3]), semantic: $(l[4]), free_ends: $(length(axon_freeends)) ($(axon_freeends))")
        if !isempty(axon_freeends) && length(axon_freeends) < 10
            #push!(checks, [l[1], axon_freeends])
            axons[l[1]] = axon_freeends
            for s in axon_freeends
                free_ends[s] = l[1]
            end
            if l[2] > 20
                long_axons[l[1]] = axon_freeends
            end
        end
        if seg_type == "axon" && l[2] > 40
            push!(really_long_axons, l[1])
        end
    elseif l[2] < 5
        if l[4][3] > l[4][1] && l[4][3] > l[4][2]
            #push!(glial_segs, l[1])
            continue
        end
        if isempty(intersect(segs[l[1]], keys(d_faceareas)))
            push!(small_pieces, l[1])
        end
    end
end
#println(glial_segs)
println("$(length(keys(axons))) axon candidates")
println("$(length(small_pieces)) small pieces")
println("$(length(keys(long_axons))) long axon candidates")
println("$(length(really_long_axons)) really long axon candidates")
#println(keys(axons))

merge_graph = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
merge_graph2 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
merge_graph3 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
considered = match_long_axons(small_pieces, long_axons, new_rg, Set{Int}(), true, merge_graph)
union!(considered, match_axons(axons, segs, new_rg, free_ends, Set{Int}(), true, merge_graph))
newly_considered = process_rg(new_rg, l_segs, segs, rg_volume, d_sizes, considered, merge_graph)
union!(considered,newly_considered)
matches = merge_edges(merge_graph, th_tier1, new_rg)

discard = match_long_axons(small_pieces, long_axons, new_rg, considered, false, merge_graph2)
union!(discard, match_axons(axons, segs, new_rg, free_ends, considered, false, merge_graph2))
matches2 = merge_edges(merge_graph2, th_tier2, new_rg)
##matches2 = merge_edges(merge_graph2, 0.199985, new_rg)
#union!(considered,discard)
#discard = match_long_axons2(long_axons, new_rg, rg_volume, segs, d_sizes, d_faceareas, considered, merge_graph3)
#matches3 = merge_edges(merge_graph3, th_tier3, new_rg)
#matches = match_long_axons(small_pieces, long_axons, new_rg)
update_new_rg(num_seg, new_rg, merge(matches, matches2))
#println([x[1] for x in checks])
