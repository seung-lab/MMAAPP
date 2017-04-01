function check_dend(segment, svInfo)
    total_vol = sum_vol(segment, svInfo)
    sem = sum_sem(segment, svInfo)
    if total_vol < 1000000 || sem[2] < sem[1] || sem[2] < 0.5*total_vol
        return false
    end
    for s in segment
        cc = check_connectivity(Set{Int}(s), segment, svInfo.regionGraph)
        if (cc > 2 && svInfo.supervoxelSizes[s] > 0.5*total_vol)
            return true
        end
    end
    return false
end

function process_edge(set_a, set_b, svInfo)
    max_aff = zero(Float64)
    rg_volume = svInfo.regionGraph
    d_sem = svInfo.semanticInfo
    for a in set_a
        for b in intersect(keys(rg_volume[a]),set_b)
            edge = rg_volume[a][b]
            if edge.num > 1000
                continue
            end
            if length(set_b) > 5 && d_sem[b][4] > 500
                continue
            end
            if edge.aff/edge.area > 0.6 || (edge.aff/edge.area > reliable_th && 20 < edge.num < 500)
                if max_aff < edge.aff/edge.area
                    println("process: $(length(set_a)), $(length(set_b))")
                    println("atomic edge: $edge")
                    max_aff = edge.aff/edge.area
                end
            end
        end
    end
    return max_aff
end

function match_segments(tails, set, trunks, rg_volume)
    max_mean = zero(Float64)
    target = zero(Int32)
    for a in tails
        for b in keys(rg_volume[a])
            if b in set || !(b in trunks)
                continue
            end
            edge = rg_volume[a][b]
            tmp = edge.aff/edge.area
            if tmp > max_mean
                max_mean = tmp
                target = b
            end
        end
    end
    if max_mean > reliable_th
        return target
    else
        return zero(Float64)
    end
end

function check_segs(segInfo, svInfo, considered, merge_graph)
    spines = Set{Int}()
    trunks = Dict{Int,Int}()
    for l in keys(segInfo.supervoxelDict)
        if check_dend(segInfo.supervoxelDict[l], svInfo)
            tails = find_ends_of_dend(segInfo.supervoxelDict[l], svInfo)
            for c in tails
                trunks[c] = l
            end
        end
    end

    for a in keys(segInfo.regionGraph)
        if a in considered
            continue
        end
        set_a = get(segInfo.supervoxelDict,a,Set{Int}(a))
        vol_a = sum_vol(set_a, svInfo)
        if length(set_a) > 100 || vol_a > 1000000
            continue
        end
        sem_a = sum_sem(set_a, svInfo)
        b, sem_max = max_sem(set_a, svInfo)
        if sem_a[4] < 200
            continue
        end
        tail = Set{Int}()
        if sem_max > 0.5*sem_a[4]
            tail = find_ends(set_a, b, svInfo.regionGraph, svInfo.semanticInfo)
        end
        current_seg = a
        while true
            seg_a, freeends_a = check_segment2(set_a, svInfo)
            ends = intersect(freeends_a, tail)
            if length(set_a) < 3
                ends = set_a
            end
            println("segid: $(current_seg), parts: $(length(set_a)), size: $(vol_a), free_ends: $(tail), $(freeends_a) ($(b))")
            if 0 < length(ends) < 3
                vol_a = sum_vol(set_a, svInfo)
                #println("segid: $(current_seg), parts: $(length(set_a)), size: $(vol_a), free_ends: $(ends) ($(b))")
                target = match_segments(ends, set_a, keys(trunks), svInfo.regionGraph)
                if target == 0 && length(set_a) > 5
                    target = match_segments(ends, set_a, keys(svInfo.regionGraph), svInfo.regionGraph)
                end
                if target != 0
                    new_seg = get(svInfo.segmentDict, target, target)
                    push!(merge_graph[current_seg], new_seg)
                    push!(merge_graph[new_seg], current_seg)
                    push!(spines, current_seg)
                    push!(spines, target)
                    println("merge: $(current_seg), $(new_seg) ($(ends), $target)")
                    new_seg_set = get(segInfo.supervoxelDict, new_seg, Set{Int}(new_seg))
                    if length(new_seg_set) < 5 && sum_vol(new_seg_set, svInfo) < 1000000
                        println("keep searching: $new_seg, $target")
                        union!(set_a, new_seg_set)
                        new_b, sem_max = max_sem(set_a, svInfo)
                        if new_b !== b
                            break
                        end
                        tail = find_ends(set_a, b, svInfo.regionGraph, svInfo.semanticInfo)
                        current_seg = new_seg
                    else
                        break
                    end
                else
                    break
                end
            else
                push!(spines, current_seg)
                break
            end
        end
    end
    println(spines)
    return spines
end

function process_rg(segInfo, svInfo, considered, merge_graph)
    visited = Set{atomic_edge}()
    dend_candidates = Set{Int}()
    spine_candidates = Set{Int}()
    trunks = Set{Int}()
    for l in keys(segInfo.supervoxelDict)
        if check_dend(segInfo.supervoxelDict[l], svInfo) && !(l in considered)
            push!(trunks, l)
        end
    end
    keep = true
    while keep
        keep = false
        for b in setdiff(keys(segInfo.regionGraph), spine_candidates)
            set_b = get(segInfo.supervoxelDict, b, Set([b]))
            if (length(set_b) > 30 || sum_vol(set_b, svInfo) > 1000000) || b in considered
                continue
            end
            sem_b = sum_sem(set_b, svInfo)
            vol_b = sum_vol(set_b, svInfo)
            if sem_b[2] < sem_b[1] || sem_b[2] < 0.3*vol_b
                continue
            end
            #if length(set_b) > 5
            #    seg_b, freeends_b = check_segment(set_b, rg_volume, d_sizes, d_facesegs)
            #    println("$b, $freeends_b")
            #    if length(freeends_b) > 0
            #        set_b = freeends_b
            #    end
            #end
            max_a = zero(Int)
            max_aff = zero(Float64)
            for a in trunks
                set_a = get(segInfo.supervoxelDict, a, Set([a]))
                #if length(set_a) > 5
                #    seg_a, freeends_a = check_segment(set_a, rg_volume, d_sizes, d_facesegs)
                #    if length(freeends_a) > 0
                #        set_a = freeends_a
                #    end
                #end
                tmp = process_edge(set_a, set_b, svInfo)
                if tmp > reliable_th && tmp > max_aff
                    max_aff = tmp
                    max_a = a
                end
            end
            if max_aff > reliable_th
                push!(dend_candidates,max_a)
                push!(spine_candidates,b)
                push!(merge_graph[max_a],b)
                push!(merge_graph[b],max_a)
                push!(segs[max_a], b)
                keep = true
            end
        end
    end
    println("dend candidates:")
    println("$dend_candidates, $(length(spine_candidates))")
    return spine_candidates
end



