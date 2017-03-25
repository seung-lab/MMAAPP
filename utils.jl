function max_sem(clst, d_sem)
    max_seg = 0
    sem_max = zero(Float32)
    for c in clst
        if d_sem[c][4] > sem_max
            sem_max = d_sem[c][4]
            max_seg = c
        end
    end
    return max_seg, sem_max
end

function sum_sem(clst, d_sem)
    total_sem = Float32[0,0,0,0]
    for c in clst
        total_sem += d_sem[c]
    end
    return total_sem
end

function max_vol(clst, d_sizes)
    max_seg = 0
    vol_max = 0
    for c in clst
        if d_sizes[c] > vol_max
            vol_max = d_sizes[c]
            max_seg = c
        end
    end
    return max_seg, vol_max
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

function merge_bbox(seg)
    bbox = [typemax(Int32),typemax(Int32),typemax(Int32),typemin(Int32),typemin(Int32),typemin(Int32)]
    for c in seg
        for i in 1:3
            bbox[i] = min(bbox[i], bboxes[c][i])
            bbox[i+3] = max(bbox[i+3], bboxes[c][i+3])
        end
    end
    return bbox
end

function expand_bbox(bbox, delta)
    for i in 1:3
        bbox[i] -= delta[i]
        bbox[i+3] += delta[i]
    end
end

function overlap_bbox(b1, b2)
    si = max(0, min(b1[4], b2[4]) - max(b1[1], b2[1])) * max(0, min(b1[5], b2[5]) - max(b1[2], b2[2])) * max(0, min(b1[6], b2[6]) - max(b1[3], b2[3]))
    if si > 0
        return true
    else
        return false
    end
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

