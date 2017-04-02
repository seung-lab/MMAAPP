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
svInfo, segInfo = load_segments([1,1,1,2048,2048,256])


segs = segInfo.supervoxelDict

axons = Dict{Int, Set{Int}}()
free_ends = Dict{Int, Int}()
small_pieces = Set{Int}()
for a in keys(segs)
    vol_a = sum_vol(segs[a], svInfo)
    sem_a = sum_sem(segs[a], svInfo)
    if vol_a < 10000000 && length(segs[a]) > 5
        #ccsz = connected(intersect(keys(rg_faces), segs[l[1]]), rg_faces, d_facesegs)
        #faces = faces_touched(segs[l[1]], face_segs)
        seg_type, axon_freeends = check_segment(segs[a], svInfo)
        if !isempty(axon_freeends) && length(axon_freeends) < 10 && (seg_type == "axon" || seg_type == "not sure")
            println("segid: $(a), parts: $(length(segs[a])), size: $(vol_a), free_ends: $(length(axon_freeends)) ($(axon_freeends))")
            #push!(checks, [l[1], axon_freeends])
            axons[a] = axon_freeends
            for s in axon_freeends
                free_ends[s] = a
            end
        end
    elseif length(segs[a]) < 5
        if isempty(intersect(segs[a], svInfo.boundarySupervoxels))
            push!(small_pieces, a)
        end
    end
end
println("$(length(keys(axons))) axon candidates")
println("$(length(small_pieces)) small pieces")
#println(keys(axons))

merge_graph1 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
merge_graph2 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
merge_graph3 = DefaultOrderedDict{Int, Set{Int}}(()->Set{Int}())
considered = Set{Int}()
union!(considered, match_long_axons(small_pieces, axons, segInfo, considered, true, merge_graph1))
union!(considered, match_axons(axons, segInfo, svInfo, free_ends, considered, true, merge_graph1))
union!(considered, process_rg(segInfo, svInfo, considered, merge_graph1))
matches = merge_edges(merge_graph1, th_tier1, segInfo.regionGraph)
union!(considered, check_segs(segInfo, svInfo, considered, merge_graph2))

#discard = match_long_axons(small_pieces, long_axons, new_rg, considered, false, merge_graph2)
#union!(discard, match_axons(axons, segs, new_rg, free_ends, considered, false, merge_graph2))
matches2 = merge_edges(merge_graph2, th_tier2, segInfo.regionGraph)
##matches2 = merge_edges(merge_graph2, 0.199985, new_rg)
#union!(considered,discard)
#discard = match_long_axons2(long_axons, new_rg, rg_volume, segs, d_sizes, d_facesegs, considered, merge_graph3)
#matches3 = merge_edges(merge_graph3, th_tier3, new_rg)
#matches = match_long_axons(small_pieces, long_axons, new_rg)
update_new_rg(svInfo.totalSupervoxels, segInfo.regionGraph, merge(matches, matches2))
#println([x[1] for x in checks])
