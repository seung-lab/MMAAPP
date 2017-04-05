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

function match_segments(tails, set, trunks, svInfo)
    max_mean = zero(Float64)
    target = zero(Int32)
    rg_volume = svInfo.regionGraph
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

function check_segs(dendrites, spines, smallSegments, segInfo, svInfo, considered, merge_graph)
    processed = Set{Int}()
    shafts = Set{Int}()
    for s in dendrites.segid
        push!(shafts, dendrites.shaft[s])
    end
    dendrite_ends = union(keys(dendrites.segment),shafts)
    for a in spines.segid
        if a in considered.segid
            continue
        end
        set_a = get(segInfo.supervoxelDict,a,Set{Int}(a))
        current_seg = a
        ends = spines.freeends[a]
        b = spines.psd[a]
        while true
            println("segid: $(current_seg), parts: $(length(set_a)), free_ends: $(ends) ($(b))")
            if 0 < length(ends) < 5
                vol_a = sum_vol(set_a, svInfo)
                target = match_segments(ends, set_a, dendrite_ends, svInfo)
                if target == 0 && length(set_a) > 5
                    target = match_segments(ends, set_a, keys(svInfo.regionGraph), svInfo)
                end
                if target != 0
                    new_seg = get(svInfo.segmentDict, target, target)
                    push!(merge_graph[current_seg], new_seg)
                    push!(merge_graph[new_seg], current_seg)
                    push!(processed, current_seg)
                    push!(processed, target)
                    println("merge: $(current_seg), $(new_seg) ($(ends), $target)")
                    new_seg_set = get(segInfo.supervoxelDict, new_seg, Set{Int}(new_seg))
                    if length(new_seg_set) < 5 && sum_vol(new_seg_set, svInfo) < 1000000
                        println("keep searching: $new_seg, $target")
                        union!(set_a, new_seg_set)
                        new_b, sem_max = max_sem(set_a, svInfo)
                        if new_b !== b
                            break
                        end
                        seg_type, freeends = check_segment(set_a, svInfo)
                        tails = find_ends(set_a, b, svInfo)
                        current_seg = new_seg
                        ends = setdiff(intersect(freeends, tails), svInfo.boundarySupervoxels)
                    else
                        break
                    end
                else
                    break
                end
            else
                push!(processed, current_seg)
                break
            end
        end
    end
    println(processed)
    return processed
end

function process_rg(axons, dendrites, smallSegments, segInfo, svInfo, considered, merge_graph)
    visited = Set{atomic_edge}()
    dend_candidates = Set{Int}()
    attached = Set{Int}()
    keep = true
    while keep
        keep = false
        for b in setdiff(smallSegments.segid, attached)
            set_b = get(segInfo.supervoxelDict, b, Set([b]))
            if b in considered.segid
                continue
            end
            sem_b = sum_sem(set_b, svInfo)
            vol_b = sum_vol(set_b, svInfo)
            if length(set_b) > 10 && !is_dendrite(sem_b, vol_b)
                continue
            end
            max_a = zero(Int)
            max_aff = zero(Float64)
            for a in dendrites.segid
                set_a = get(segInfo.supervoxelDict, a, Set{Int}(a))
                tmp = process_edge(set_a, set_b, svInfo)
                if tmp > reliable_th && tmp > max_aff
                    max_aff = tmp
                    max_a = a
                end
            end
            if max_aff > reliable_th
                push!(dend_candidates,max_a)
                push!(attached,b)
                push!(merge_graph[max_a],b)
                push!(merge_graph[b],max_a)
                println("merge: $(b), $(max_a)")
                set_a = get(segInfo.supervoxelDict, max_a, Set{Int}(max_a))
                push!(set_a, b)
                segInfo.supervoxelDict[max_a] = set_a
                keep = true
            end
        end
    end
    println("dend candidates:")
    println("$dend_candidates, $(length(attached))")
    return attached
end



