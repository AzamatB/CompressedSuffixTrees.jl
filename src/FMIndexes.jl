module FMIndexes

export FMIndex, restore, count, locate, locateall

include("WaveletMatrices.jl")

using IndexableBitVectors
using Mmap
using SuffixArrays

using .WaveletMatrices

"""
Index for full-text search.

Type parameters:

* `W`: the number of bits required to encode the alphabet
* `T`: the type to represent positions of a sequence
"""
struct FMIndex{W,T}
    bwt::WaveletMatrix{W,UInt16,SucVector}
    sentinel::Int
    suffixes::Vector{T}
    count::Vector{Int}
end

function FMIndex(seq, suffixes, σ::Integer)
    wm = WaveletMatrix(make_bwt(seq, suffixes), log2(Int, σ))
    sentinel = something(findfirst(isequal(0), suffixes), 0) + 1
    # count characters
    count = count_bytes(seq, σ)
    count[1] = 1  # sentinel '$' is smaller than any character
    cumsum!(count, count)
    return FMIndex(wm, sentinel, suffixes, count)
end

"""
    FMIndex(seq, σ)

Build an FM-Index from a sequence `seq`.
The sequence must support `convert(UInt, seq[i])` for each character and the
alphabet size should be less than or equal to 256. The second parameter, `σ`, is
the alphabet size.
"""
function FMIndex(seq, σ::Integer=(maximum(seq) + 2))
    T = index_type(length(seq))
    @assert 1 ≤ σ ≤ typemax(UInt16) + 1
    suffixes = make_sa(T, seq, σ)
    return FMIndex(seq, suffixes, σ)
end

"""
    FMIndex(text)

Build an FM-Index from an ASCII text.
"""
function FMIndex(text::Union{String,SubString{String}})
    return FMIndex(codeunits(text), 128)
end

Base.length(index::FMIndex) = length(index.bwt)

"""
Restore the original text from the index.
"""
function restore(index::FMIndex)
    n = length(index)
    text = Vector{UInt16}(undef, n)
    p = index.sentinel
    while n > 0
        p = lfmap(index, p)
        text[n] = index.bwt[p ≥ index.sentinel ? p - 1 : p]
        n -= 1
    end
    return text
end

"""
Count the number of occurrences of the given query.
"""
function Base.count(query, index::FMIndex)
    return length(sa_range(query, index))
end

"""
Locate the positions of occurrences of the given query.
This method returns an iterator of positions:

    for pos in locate(query, index)
        # ...
    end
"""
function locate(query, index::FMIndex)
    return LocationIterator(sa_range(query, index), index)
end

"""
Locate the positions of all occurrences of the given query.
"""
function locateall(query, index::FMIndex)
    iter = locate(query, index)
    locs = Vector{Int}(undef, length(iter))
    for (i, loc) in enumerate(iter)
        locs[i] = loc
    end
    return locs
end

function log2(::Type{Int}, x::Integer)
    return 64 - leading_zeros(convert(UInt64, x - 1))
end

function count_bytes(seq, σ)
    count = zeros(Int, σ + 1)
    for i in 1:length(seq)
        count[convert(UInt16, seq[i])+2] += 1
    end
    resize!(count, σ)
    return count
end

# LF-mapping
function lfmap(index::FMIndex, i)
    if i == index.sentinel
        return 1
    elseif i > index.sentinel
        i -= 1
    end
    char = index.bwt[i]
    @inbounds return index.count[char+1] + rank(char, index.bwt, i)
end

function sa_range(query, index::FMIndex)
    sa_range(query, index::FMIndex, 1:(length(index)+1))
end

function sa_range(query, index::FMIndex, init_range::UnitRange{Int})
    sp, ep = init_range.start, init_range.stop
    # backward search
    i = length(query)
    while sp ≤ ep && i ≥ 1
        char = convert(UInt, query[i])
        c = index.count[char+1]
        sp = c + rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
        i -= 1
    end
    return sp:ep
end

@inline function sa_range(char::UInt, index::FMIndex, range::UnitRange{Int})
    sp, ep = range.start, range.stop
    c = index.count[char+1]
    sp = c + rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
    ep = c + rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
    return sp:ep
end

function sa_value(i::Int, index::FMIndex)
    if i == 1
        # point to the sentinel '$'
        return length(index) + 1
    end
    return index.suffixes[i - 1]
end

sa_value(i::Integer, index::FMIndex) = sa_value(Int(i), index)

# Suffix Array Construction Algorithms

# a wrapper type for a sequence that returns byte-convertible elements
struct ByteSeq{S} <: AbstractVector{UInt}
    data::S
end
Base.size(seq::ByteSeq) = (length(seq.data),)
@inline Base.getindex(seq::ByteSeq, i::Integer) = UInt(seq.data[i])

# SuffixArrays.jl: https://github.com/quinnj/SuffixArrays.jl
function make_sa(T, seq, σ)
    n = length(seq)
    suffixes_temp = Vector{Int}(undef, n)
    SuffixArrays.sais(ByteSeq(seq), suffixes_temp, 0, n, nextpow(2, σ), false)
    suffixes = Vector{T}(undef, n)
    copyto!(suffixes, suffixes_temp)
    return suffixes
end

# other utils
function index_type(n)
    n -= 1
    n ≤ typemax(UInt8) ? UInt8 :
    n ≤ typemax(UInt16) ? UInt16 :
    n ≤ typemax(UInt32) ? UInt32 : UInt64
end

# Burrows-Wheeler Transform
function make_bwt(seq, sa)
    n = length(seq)
    @assert length(sa) == n
    ret = Vector{UInt16}(undef, n)
    j = 1
    for i in 1:n
        # note that `sa` starts from zero
        p = sa[i]
        if p == 0
            ret[1] = seq[end]
        else
            ret[j+=1] = seq[p]
        end
    end
    return ret
end

struct LocationIterator{W}
    range::UnitRange{Int}
    index::FMIndex{W}
end

Base.length(iter::LocationIterator) = length(iter.range)

function Base.iterate(iter::LocationIterator, i::Int=1)
    if i > length(iter)
        return nothing
    end
    return sa_value(iter.range[i], iter.index) + 1, i + 1
end

end # FMIndexes
