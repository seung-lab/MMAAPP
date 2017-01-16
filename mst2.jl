using HDF5
using EMIRT
using DataStructures

type MST
    dend::Array{Array{UInt32,1},1}
    dendValues::Array{Float32,1}
end

function processln(l)
    data = split(l)
    parent = parse(Int, data[3])
    child1 = parse(Int, data[1])
    child2 = parse(Int, data[2])
    weight = parse(Float32, data[4])
    u1 = parse(Int, data[5])
    u2 = parse(Int, data[6])
    if parent == child1
        return (parent, child2, weight, u1, u2)
    else
        return (parent, child1, weight, u1, u2)
    end
end

function graph_to_tree(region_graph)
    #BFS
    tree = MST(Array{UInt32,1}[], Array{Float32,1}())
    visited = Set()
    roots = keys(region_graph)
    println("len: $(length(roots))")
    count = 0
    for root in roots
        if root in visited
            continue
        end
        count += 1
        queue = Queue(Int)
        enqueue!(queue, root)

        while length(queue) > 0
            root = dequeue!(queue)
            if root in visited
                continue
            end

            push!(visited,root)
            for (neighboor, weight) in region_graph[root]
                if !(neighboor in visited)
                    push!(tree.dend, UInt32[UInt32(root), UInt32(neighboor)])
                    push!(tree.dendValues, Float32(weight))
                    enqueue!(queue, neighboor)
                end
            end
        end
    end
    println(count)

    return tree

end

rootNodes = Dict{UInt32, Float32}()

mst = MST(Array{UInt32,1}[], Array{Float32,1}())


atomic_region_graph = DefaultDict(Int , Dict{Int , Real},
                   ()->Dict{Int , Real}())

open(ARGS[1]) do f
    for ln in eachline(f)
        parent, child, weight, head, tail = processln(ln)
        #println("$parent, $child, $weight")
        min_weight = weight
        if haskey(rootNodes, child)
            min_weight = min(min_weight, rootNodes[child])
            delete!(rootNodes,child)
        end
        if haskey(rootNodes, parent)
            min_weight = min(min_weight, rootNodes[parent])
        end
        rootNodes[parent] = min_weight
        atomic_region_graph[head][tail] = min_weight
        atomic_region_graph[tail][head] = min_weight
    end
end

mst = graph_to_tree(atomic_region_graph)

N=length(mst.dendValues)
println("length: $N")
dend = zeros(UInt32, N,2)
for i in 1:N
    dend[i,:]=mst.dend[i]
end
dend[:,1], dend[:,2] = dend[:,2], dend[:,1]
f = h5open(ARGS[2], "r+")
f["dend"] = dend
f["dendValues"] = mst.dendValues
close(f)
