using Base.Threads
using DataStructures
using Combinatorics
include("segments.jl")
include("output.jl")
include("utils.jl")
include("classify.jl")
include("freeends.jl")
include("spine.jl")
include("axon.jl")

agg_threshold = parse(Float64, ARGS[1])
th_tier1 = agg_threshold - 0.000005
th_tier2 = agg_threshold - 0.000015
th_tier3 = agg_threshold - 0.000025
reliable_th = parse(Float64, ARGS[2])

svInfo = SupervoxelInfo()
segInfo = SegmentInfo()

@time load_segments!([28673, 24065, 2, 80634, 60859, 996], svInfo, segInfo)

axons = Axons()
dendrites = Dendrites()
spines = Spines()
smallSegments = SmallSegments()
processedSegments = ProcessedSegments()

@time classify_segments!(segInfo, svInfo, axons, dendrites, spines, smallSegments, processedSegments)


merge_graph1 = DefaultOrderedDict{UInt64, Set{UInt64}}(()->Set{UInt64}())
merge_graph2 = DefaultOrderedDict{UInt64, Set{UInt64}}(()->Set{UInt64}())
merge_graph3 = DefaultOrderedDict{UInt64, Set{UInt64}}(()->Set{UInt64}())
considered = Set{UInt64}()
#union!(considered, match_long_axons(small_pieces, axons, segInfo, considered, true, merge_graph1))
@time union!(processedSegments.segid, match_axons(axons, smallSegments, segInfo, svInfo, processedSegments, merge_graph1))
@time union!(processedSegments.segid, process_rg(axons, dendrites, smallSegments, segInfo, svInfo, processedSegments, merge_graph1))
matches = merge_edges(merge_graph1, th_tier1, segInfo.regionGraph)
@time union!(processedSegments.segid, check_segs(dendrites, spines, smallSegments, segInfo, svInfo, processedSegments, merge_graph2))

#discard = match_long_axons(small_pieces, long_axons, new_rg, considered, false, merge_graph2)
#union!(discard, match_axons(axons, segs, new_rg, free_ends, considered, false, merge_graph2))
matches2 = merge_edges(merge_graph2, th_tier2, segInfo.regionGraph)
##matches2 = merge_edges(merge_graph2, 0.199985, new_rg)
#union!(considered,discard)
#discard = match_long_axons2(long_axons, new_rg, rg_volume, segs, d_sizes, d_facesegs, considered, merge_graph3)
#matches3 = merge_edges(merge_graph3, th_tier3, new_rg)
#matches = match_long_axons(small_pieces, long_axons, new_rg)
@time update_new_rg(svInfo.totalSupervoxels, segInfo.regionGraph, merge(matches, matches2))
#println([x[1] for x in checks])
