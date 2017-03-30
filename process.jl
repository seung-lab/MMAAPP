using HDF5
using EMIRT
using DataStructures
using Combinatorics
include("segments.jl")
include("output.jl")
include("utils.jl")
include("freeends.jl")
include("spine.jl")
include("axon.jl")

agg_threshold = parse(Float64, ARGS[1])
th_tier1 = agg_threshold - 0.000005
th_tier2 = agg_threshold - 0.000015
th_tier3 = agg_threshold - 0.000025
reliable_th = parse(Float64, ARGS[2])

segs, pd = agglomerate("test_mst.in")
num_seg, rg_volume = read_rg("rg_volume.in", Dict{UInt32,UInt32}())
d_sizes = read_size("sv_size.in")
d_sem = read_semantic("sem_volume.in")
#rg_faces, face_segs, d_facesegs = process_faces(sgm.segmentation)
_, new_rg = read_rg("new_rg.in", pd)
bboxes = read_bboxes("bbox_volume.in")
println("$(length(keys(d_sizes)))")
println("size of rg: $(length(keys(new_rg)))")
#d_facesegs = DefaultOrderedDict{Int,Int}(0)
d_facesegs = find_facesegs(bboxes,[1,1,1,2048,2048,256])

l_segs = []

for k in keys(segs)
    push!(l_segs, [k,length(segs[k]),sum_vol(segs[k], d_sizes),sum_sem(segs[k], d_sem)])
end

#sort!(l_segs, by=x->x[3]/x[2],rev=true)

#println("size of rg: $(length(keys(rg_faces)))")
axons = Dict{Int, Set{Int}}()
free_ends = Dict{Int, Int}()
small_pieces = Set{Int}()
long_axons = Dict{Int, Set{Int}}()
really_long_axons = Set{Int}()
for l in l_segs
    if l[3] < 10000000 && l[2] > 5
        #ccsz = connected(intersect(keys(rg_faces), segs[l[1]]), rg_faces, d_facesegs)
        #faces = faces_touched(segs[l[1]], face_segs)
        seg_type, axon_freeends = check_segment(segs[l[1]], rg_volume, d_sizes, d_facesegs)
        if !isempty(axon_freeends) && length(axon_freeends) < 10
            println("segid: $(l[1]), parts: $(l[2]), size: $(l[3]), free_ends: $(length(axon_freeends)) ($(axon_freeends))")
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
        if isempty(intersect(segs[l[1]], d_facesegs))
            push!(small_pieces, l[1])
        end
    end
end
println("$(length(keys(axons))) axon candidates")
println("$(length(small_pieces)) small pieces")
println("$(length(keys(long_axons))) long axon candidates")
println("$(length(really_long_axons)) really long axon candidates")
#println(keys(axons))

merge_graph1 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
merge_graph2 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
merge_graph3 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
considered = Set{Int}()
union!(considered, match_long_axons(small_pieces, long_axons, new_rg, considered, true, merge_graph1))
union!(considered, match_axons(axons, segs, new_rg, free_ends, considered, true, merge_graph1))
union!(considered, process_rg(new_rg, segs, rg_volume, d_sizes, d_sem, considered, merge_graph1))
matches = merge_edges(merge_graph1, th_tier1, new_rg)
union!(considered, check_segs(new_rg, segs, rg_volume, d_sizes, d_sem, considered, merge_graph2))

#discard = match_long_axons(small_pieces, long_axons, new_rg, considered, false, merge_graph2)
#union!(discard, match_axons(axons, segs, new_rg, free_ends, considered, false, merge_graph2))
matches2 = merge_edges(merge_graph2, th_tier2, new_rg)
##matches2 = merge_edges(merge_graph2, 0.199985, new_rg)
#union!(considered,discard)
#discard = match_long_axons2(long_axons, new_rg, rg_volume, segs, d_sizes, d_facesegs, considered, merge_graph3)
#matches3 = merge_edges(merge_graph3, th_tier3, new_rg)
#matches = match_long_axons(small_pieces, long_axons, new_rg)
update_new_rg(num_seg, new_rg, merge(matches, matches2))
#println([x[1] for x in checks])
