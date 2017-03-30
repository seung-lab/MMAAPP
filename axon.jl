function check_edge(rg_volume, edge, seg1, seg2)
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

function match_axons(axons, segs, new_rg, rg_volume, free_ends, considered, is_strict, merge_graph)
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
            if check_edge(rg_volume, a_edge, segs[a], segs[b])
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

function match_long_axons2(long_axons, new_rg, rg_volume, segs, d_sizes, d_facesegs, considered, merge_graph)
    pairs = []
    processed = Set{Int}()
    for s in keys(d_sizes)
        if s in d_facesegs
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
