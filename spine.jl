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
    println("targets: $(keys(dendrites.segment))")
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
                target = match_segments(ends, set_a, keys(dendrites.segment), svInfo)
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



