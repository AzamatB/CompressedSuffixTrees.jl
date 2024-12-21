const ∞ = typemax(Int)

mutable struct Node{T}
    children::Dict{T,Int}
    start::Int
    ending::Int
    suffix_link::Int
    suffix_index::Int

    function Node{T}(start::Integer=0, ending::Integer=∞) where {T}
        return new{T}(Dict{T,Int}(), start, ending, 0, -1)
    end
end


mutable struct CompressedSuffixTree{T}
    nodes::Vector{Node{T}}
    text::Vector{T}
    root::Int
    position::Int
    current_node::Int
    suffix_link_node::Int
    remainder::Int
    active_node::Int
    active_len::Int
    active_edge::Int
end

function add_suffix_link(cst::CompressedSuffixTree, node_idx::Integer)
    if cst.suffix_link_node > 0
        node = cst.nodes[cst.suffix_link_node]
        node.suffix_link = node_idx
    end
    cst.suffix_link_node = node_idx
    return node_idx
end

function edge_length(cst::CompressedSuffixTree, node::Node)
    len = min(node.ending, cst.position + 1) - node.start
    return len
end

function add_node!(cst::CompressedSuffixTree, start::Integer, ending::Integer)
    node_idx = cst.current_node + 1
    cst.current_node = node_idx
    cst.nodes[node_idx] = Node(start, ending)
    return node_idx
end



Node{Char}()
