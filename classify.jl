type Axons
    segid::SupervoxelSet
    segment::SegmentDict # freeend->segid
    freeends::SupervoxelDict # segid->freeends
    function Axons()
        s1 = SupervoxelSet()
        s2 = SegmentDict()
        f = SupervoxelDict()
        return new(s1,s2,f)
    end
end

type Dendrites
    segid::SupervoxelSet
    shaft::SegmentDict # segid->shaft
    segment::SegmentDict # freeend->segid
    freeends::SupervoxelDict # segid->freeends
    anchor::SegmentDict # freeend->anchor
    function Dendrites()
        s1 = SupervoxelSet()
        s2 = SegmentDict()
        s3 = SegmentDict()
        f = SupervoxelDict()
        a = SegmentDict()
        return new(s1,s2,s3,f,a)
    end
end

type Spines
    segid::SupervoxelSet
    psd::SegmentDict # segid->psd
    segment::SegmentDict # freeend->segid
    freeends::SupervoxelDict # segid->freeends
    anchor::SegmentDict # freeend->anchor
    function Spines()
        s1 = SupervoxelSet()
        p = SegmentDict()
        s2 = SegmentDict()
        f = SupervoxelDict()
        a = SegmentDict()
        return new(s1,p,s2,f,a)
    end
end

type SmallSegments
    segid::SupervoxelSet
    function SmallSegments()
        return new(SupervoxelSet())
    end
end

type ProcessedSegments
    segid::SupervoxelSet
    function ProcessedSegments()
        return new(SupervoxelSet())
    end
end

const vol_threshold = 1000000

function is_glial(sem, vol)
    return sem[3] > 0.5*vol
end

function is_dendrite(sem, vol)
    return sem[2] > 0.5*vol && sem[2] > sem[1]
end

function is_axon(sem, vol)
    return sem[1] > 0.5*vol && sem[1] > sem[2]
end

function find_psd(seg, svInfo)
    sem = sum_sem(seg, svInfo)
    b, sem_max = max_sem(seg, svInfo)
    if sem[4] < 200
        return 0
    end
    if sem_max > 0.5*sem[4]
        return b
    end
    return 0
end

function check_semantic(segment, svInfo)
    sem = sum_sem(segment, svInfo)
    vol = sum_vol(segment, svInfo)
    seg_type = "not sure"
    if is_glial(sem, vol)
        seg_type = "glial"
    elseif is_dendrite(sem, vol)
        seg_type = "dendrite"
    elseif is_axon(sem, vol)
        seg_type = "axon"
    else
        seg_type = "not sure"
    end
    return seg_type
end


function check_segment(segment, svInfo)
    free_ends = Set{UInt64}()
    seg_type = ""
    if length(segment) < 5
        return seg_type, free_ends
    end
    total_vol = sum_vol(segment, svInfo)
    largest_s, vol_max = max_vol(segment, svInfo)
    seg_type = check_semantic(segment, svInfo)
    cc = check_connectivity(Set{UInt64}(largest_s), segment, svInfo.regionGraph)
    if total_vol > vol_threshold && (cc > 2 && vol_max > 0.5*total_vol)
        if seg_type == "not sure"
            seg_type = "dendrite"
        end
        return seg_type, free_ends
    end
    #seg_type = "not sure"
    free_ends = find_ends(segment, largest_s, svInfo)
    if cc == 1
        push!(free_ends, largest_s)
    end
    #ends = setdiff(free_ends, keys(d_facesegs))
    ends = really_check_freeends(free_ends, segment, svInfo)
    #if seg_type == "not sure"
    #    if length(ends) > 10
    #        if count * 3 < length(segment)
    #            seg_type = "glial"
    #        else
    #            seg_type = "dendrite"
    #        end
    #    elseif count * 2 > length(segment)
    #        seg_type = "axon"
    #    end
    #end
    return seg_type, ends
end


function classify_segments(segInfo, svInfo)
    axons = Axons()
    dendrites = Dendrites()
    spines = Spines()
    smallSegments = SmallSegments()
    processedSegments = ProcessedSegments()

    segs = segInfo.supervoxelDict
    segids = collect(keys(segs))
    num_seg = length(segids)
    seg_types = Array{String}(num_seg)
    freeendss = Array{Set{UInt64}}(num_seg)
    @threads for i in 1:num_seg
        seg_types[i], freeendss[i] = check_segment(segs[segids[i]], svInfo)
    end
    seg_type_dict = Dict{UInt64, String}()
    freeends_dict = Dict{UInt64, Set{UInt64}}()
    for i in 1:num_seg
        seg_type_dict[segids[i]] = seg_types[i]
        freeends_dict[segids[i]] = freeendss[i]
    end
    for a in keys(segs)
        vol_a = sum_vol(segs[a], svInfo)
        sem_a = sum_sem(segs[a], svInfo)
        if length(segs[a]) >= 5
            seg_type, freeends = check_segment(segs[a], svInfo)
            println("segid: $(a), parts: $(length(segs[a])), size: $(vol_a), free_ends: $(length(freeends_dict[a])) ($(freeends_dict[a])) $(seg_type_dict[a])")
            if seg_type_dict[a] == "glial" && length(segs[a]) > 30
                push!(processedSegments.segid, a)
                continue
            end
            if length(segs[a]) <= 30 && vol_a < vol_threshold
                push!(smallSegments.segid, a)
            end
            if !isempty(freeends_dict[a]) && (seg_type_dict[a] == "axon" || seg_type_dict[a] == "not sure")
                println("axon: segid: $(a), parts: $(length(segs[a])), size: $(vol_a), free_ends: $(length(freeends_dict[a])) ($(freeends_dict[a]))")
                push!(axons.segid, a)
                axons.freeends[a] = freeends_dict[a]
                for s in freeends_dict[a]
                    axons.segment[s] = a
                end
            end
            if vol_a > vol_threshold && seg_type_dict[a] == "dendrite"
                seg = segs[a]
                m, _ = max_vol(seg, svInfo)
                push!(dendrites.segid, a)
                dendrites.shaft[a] = m
                branches = setdiff(seg, Set{UInt64}(m))
                anchor = intersect(branches, keys(svInfo.regionGraph[m]))
                freeends = SupervoxelSet()
                anchor_array = collect(anchor)
                tails_array = Array{Set{UInt64}}(length(anchor_array))
                for i in 1:length(anchor_array)
                    tails_array[i] = setdiff(find_ends(branches,anchor_array[i],svInfo), svInfo.boundarySupervoxels)
                end
                for i in 1:length(anchor_array)
                    union!(freeends, tails_array[i])
                    for t in tails_array[i]
                        dendrites.segment[t] = a
                        dendrites.anchor[t] = anchor_array[i]
                    end
                end
                println("dendrite: segid: $(a), parts: $(length(segs[a])), size: $(vol_a), free_ends: $(length(freeends)) ($(freeends)) $(seg_type_dict[a])")
                dendrites.freeends[a] = freeends
            elseif length(segs[a]) < 100 && vol_a < vol_threshold
                b = find_psd(segs[a], svInfo)
                if b == 0
                    continue
                end
                println("spine: segid: $(a), parts: $(length(segs[a])), size: $(vol_a), free_ends: $(length(freeends_dict[a])) ($(freeends_dict[a])) $(seg_type_dict[a])")
                push!(spines.segid, a)
                spines.psd[a] = b
                tails = find_ends(segs[a], b, svInfo;free=false)
                spines.freeends[a] = setdiff(intersect(freeends_dict[a], tails), svInfo.boundarySupervoxels)
                for t in spines.freeends[a]
                    spines.segment[t] = a
                    spines.anchor[t] = b
                end
            end
        elseif length(segs[a]) < 5 && vol_a < vol_threshold
            #if is_glial(sum_sem(segs[a], svInfo), sum_vol(segs[a], svInfo))
            #    push!(processedSegments.segid, a)
            #    continue
            #end
            push!(smallSegments.segid, a)
            b = find_psd(segs[a], svInfo)
            if b != 0
                push!(spines.segid, a)
                spines.psd[a] = b
                spines.freeends[a] = segs[a]
                for t in spines.freeends[a]
                    spines.segment[t] = a
                    spines.anchor[t] = b
                end
            end
        end
    end
    for a in setdiff(keys(segInfo.regionGraph), keys(segs))
        #if is_glial(svInfo.semanticInfo[a], svInfo.supervoxelSizes[a])
        #    push!(processedSegments.segid, a)
        #    continue
        #end
        push!(smallSegments.segid, a)
        if svInfo.semanticInfo[a][4] > 200
            push!(spines.segid, a)
            spines.psd[a] = a
            spines.freeends[a] = SupervoxelSet(a)
            spines.anchor[a] = a
        end
    end
    seg_types = []
    freeendss = []
    gc()

    return axons, dendrites, spines, smallSegments, processedSegments
end

