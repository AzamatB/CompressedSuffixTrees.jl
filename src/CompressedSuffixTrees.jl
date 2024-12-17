module CompressedSuffixTrees

include("FMIndexes.jl")

using .FMIndexes

"""
Kasai's algorithm for constructing the LCP array.
Takes a string and its suffix array as input.
Returns the LCP array where LCP[i] is the length of the longest common prefix
between suffixes starting at SA[i] and SA[i-1].
Time complexity: O(n)
"""
function build_lcp_array(text::Vector{T}, suffix_array::Vector{Int}) where {T<:Unsigned}
    n = length(text)
    lcp = zeros(Int, n)
    rank = zeros(Int, n)

    # Compute rank array (inverse of suffix array)
    for i in 1:n
        rank[suffix_array[i]] = i
    end
    # Initialize LCP length
    k = 0
    # Compute LCP values
    for i in 1:n
        if rank[i] == n
            k = 0
            continue
        end
        j = suffix_array[rank[i]+1]
        # Extend the previous LCP value
        while (i + k ≤ n) && (j + k ≤ n) && (text[i+k] == text[j+k])
            k += 1
        end
        lcp[rank[i]] = k
        # Update k for the next iteration
        k = max(0, k - 1)
    end
    return lcp
end


struct CompressedSuffixTree{T}
    text::Vector{T}
    fm_index::FMIndex{T}
end

function CompressedSuffixTree(text::AbstractVector)
    fm_index = FMIndex(text)
    return CompressedSuffixTree(text, fm_index)
end

# Helper function to get the range of suffixes prefixed by a pattern
function suffix_range(cst::CompressedSuffixTree{T}, pattern::AbstractVector{T}) where {T}
    range = fm_index_range(cst.fm_index, pattern)
    return range
end

# Function to check if a pattern exists
function contains(cst::CompressedSuffixTree{T}, pattern::Vector{T}) where {T}
    range = suffix_range(cst, pattern)
    return !isempty(range)
end


end # CompressedSuffixTrees
