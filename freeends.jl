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

function really_check_freeends(ends, segment, svInfo)
    free_ends = Set{Int}()
    boundary_segs = svInfo.boundarySupervoxels
    rg_volume = svInfo.regionGraph
    avg_vol = sum_vol(segment, svInfo)/length(segment)
    for s in ends
        end_segment = Set{Int}()
        near_boundary = false
        if s in boundary_segs
            continue
        end
        #if d_sizes[s] > 5*avg_vol
        #    continue
        #end
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

function check_segment(segment, svInfo)
    free_ends = Set{Int}()
    count = 0
    total_vol = sum_vol(segment, svInfo)
    seg_type = "not sure"
    for s in segment
        cc = check_connectivity(Set{Int}(s), segment, svInfo.regionGraph)
        if total_vol > 1000000 && (cc > 2 && svInfo.supervoxelSizes[s] > 0.5*total_vol)
            seg_type = "dend"
        end
        if cc == 2
            count += 1
        elseif cc == 1
            push!(free_ends, s)
        end
    end
    #ends = setdiff(free_ends, keys(d_facesegs))
    ends = really_check_freeends(free_ends, segment, svInfo)
    if seg_type == "not sure" && length(ends) > 3 && (count * 3 < length(segment))
        seg_type = "glial"
    end
    if count * 2 > length(segment)
        return "axon", ends
    else
        return seg_type, ends
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

function find_ends(segs, head, rg_volume, d_sem)
    visited = Set{Int}()
    queue = Queue(Int)
    enqueue!(queue, head)
    tail = Set{Int}([head])
    children = Set{Int}()
    while length(queue) > 0
        root = dequeue!(queue)
        if !(root in visited)
            if root != head && d_sem[root][4] > 500
                return Set{Int}()
            end
            push!(visited, root)
            for neighboor in keys(rg_volume[root])
                if (!(neighboor in visited)) && neighboor in segs
                    push!(children, neighboor)
                end
            end
        end
        if length(queue) == 0 && length(children) > 0
            for c in children
                enqueue!(queue, c)
            end
            tail = children
            children = Set{Int}()
        end
    end
    return tail
end

function find_ends_of_dend(seg, svInfo)
    vol_seg = sum_vol(seg, svInfo)
    a, vol_max = max_vol(seg, svInfo)
    if vol_max < 0.5*vol_seg
        return Set{Int}()
    end
    branches = setdiff(seg, Set{Int}(a))
    tails = Set{Int}()
    for b in intersect(branches, keys(svInfo.regionGraph[a]))
        union!(tails, find_ends(branches,b,svInfo.regionGraph,svInfo.semanticInfo))
    end
    return tails
end


