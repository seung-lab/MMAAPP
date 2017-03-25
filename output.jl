function print_edge(rg_out, edge, mean_aff)
    area = edge.num
    sum_aff = edge.sum
    if mean_aff > 0
        sum_aff = area*mean_aff
    else
        sum_aff *= 0.9995
    end
    p = minmax(edge.p1, edge.p2)
    push!(rg_out,"$(p[1]) $(p[2]) $(sum_aff) $(area) $(edge.v1) $(edge.v2) $(edge.aff) $(edge.area)\n")
end

function update_new_rg(num_seg, new_rg, matches)
    f = open("axon.in","w")
    visited = Set{atomic_edge}()
    rg_out = AbstractString[]
    for edge in keys(matches)
        print_edge(rg_out, edge, matches[edge])
        push!(visited, edge)
    end
    for a in keys(new_rg)
        for b in keys(new_rg[a])
            edge = new_rg[a][b]
            if edge in visited
                continue
            else
                print_edge(rg_out, edge, -1)
                push!(visited, edge)
            end
        end
    end
    write(f,"$(num_seg-1) $(num_seg) $(length(rg_out))\n")
    for str in rg_out
        write(f, str)
    end
    close(f)
end

function merge_edges(merge_graph, threshold, new_rg)
    edges = OrderedDict{atomic_edge, Float64}()
    visited = Set{Int}()
    for p in keys(merge_graph)
        axons_group = Int[]
        if p in visited
            continue
        end
        queue = Queue(Int)
        enqueue!(queue, p)
        push!(axons_group, p)

        while length(queue) > 0
            root = dequeue!(queue)
            if root in visited
                continue
            end

            push!(visited,root)
            for neighboor in merge_graph[root]
                if !(neighboor in visited)
                    enqueue!(queue, neighboor)
                    push!(axons_group, neighboor)
                end
            end
        end
        for q in combinations(axons_group,2)
            if haskey(new_rg[q[1]], q[2])
                edges[new_rg[q[1]][q[2]]] = threshold
            end
        end
    end
    return edges
end

