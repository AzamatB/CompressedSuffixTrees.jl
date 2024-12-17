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
    bwt::WaveletMatrix{W,UInt,SucVector}
    sentinel::Int
    samples::Vector{T}
    sampled::SucVector
    count::Vector{Int}
end

function FMIndex(seq, sa, σ, r)
    wm = WaveletMatrix(make_bwt(seq, sa), log2(Int, σ))
    # sample suffix array
    samples, sampled = sample_sa(sa, r)
    sentinel = something(findfirst(isequal(0), sa), 0) + 1
    # count characters
    count = count_bytes(seq, σ)
    count[1] = 1  # sentinel '$' is smaller than any character
    cumsum!(count, count)
    return FMIndex(wm, sentinel, samples, SucVector(sampled), count)
end

"""
    FMIndex(seq, σ=256; r=32, program=:SuffixArrays, mmap::Bool=false, opts...)

Build an FM-Index from a sequence `seq`.
The sequence must support `convert(UInt, seq[i])` for each character and the
alphabet size should be less than or equal to 256. The second parameter, `σ`, is
the alphabet size. The third parameter, `r`, is the interval of sampling values
from a suffix array. If you set it large, you can save the memory footprint but
it requires more time to locate the position.
"""
function FMIndex(seq, σ=256; r=32, program=:SuffixArrays, mmap::Bool=false, opts...)
    T = index_type(length(seq))
    opts = Dict(opts)
    @assert !haskey(opts, :σ) "σ should be passed as the second argument"
    if program === :SuffixArrays
        @assert 1 ≤ σ ≤ typemax(UInt) + 1
        sa = make_sa(T, seq, σ, mmap)
    elseif program === :psascan || program === :pSAscan
        @assert 1 ≤ σ ≤ typemax(UInt)
        psascan = get(opts, :psascan, "psascan")
        workdir = get(opts, :workdir, pwd())
        sa = make_sa_pscan(T, seq, psascan, workdir, mmap)
    else
        error("unknown program name: $program")
    end
    return FMIndex(seq, sa, σ, r)
end

"""
    FMIndex(text; opts...)

Build an FM-Index from an ASCII text.
"""
function FMIndex(text::Union{String,SubString{String}}; opts...)
    return FMIndex(codeunits(text), 128; opts...)
end

Base.length(index::FMIndex) = length(index.bwt)

function Base.show(io::IO, fmindex::FMIndex)
    println(io, summary(fmindex), ':')
    totalsize = (
        sizeof(fmindex.bwt) +
        sizeof(fmindex.samples) +
        sizeof(fmindex.sampled) +
        sizeof(fmindex.count)
    )
    print("     length: ", length(fmindex), '\n')
    print("  data size: ", Humanize.datasize(totalsize, style=:bin))
end

"""
Restore the original text from the index.
"""
function restore(index::FMIndex)
    n = length(index)
    text = Vector{UInt}(undef, n)
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
        count[convert(UInt, seq[i])+2] += 1
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
    d = 0
    @inbounds while !index.sampled[i-1]
        i = lfmap(index, i)
        d += 1
    end
    return index.samples[rank1(index.sampled, i - 1)] + d
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
function make_sa(T, seq, σ, mmap)
    n = length(seq)
    tmp_sa = mmap ? Mmap.mmap(Vector{Int}, n) : Vector{Int}(undef, n)
    SuffixArrays.sais(ByteSeq(seq), tmp_sa, 0, n, nextpow(2, σ), false)
    sa = mmap ? Mmap.mmap(Vector{T}, n) : Vector{T}(undef, n)
    copyto!(sa, tmp_sa)
    return sa
end

# pSAscan: https://www.cs.helsinki.fi/group/pads/pSAscan.html
function make_sa_pscan(T, seq, psascan, workdir, mmap)
    seqpath, io = mktemp(workdir)
    sapath = string(seqpath, ".sa5")
    try
        dump_seq(io, seq)
        run(`$psascan -o $sapath $seqpath`)
        return load_sa(T, sapath, mmap)
    catch
        rethrow()
    finally
        rm(seqpath)
        isfile(sapath) && rm(sapath)
    end
end

function dump_seq(io, seq)
    @inbounds for i in 1:length(seq)
        write(io, convert(UInt, seq[i]))
    end
    close(io)
end

function load_sa(T, file, mmap)
    # load a 40-bit suffix array generated from psascan
    size = filesize(file)
    @assert size % 5 == 0 "file $file is not 40-bit integers"
    n = div(size, 5)
    sa = mmap ? Mmap.mmap(Vector{T}, n) : Vector{T}(n)
    open(file) do input
        load_sa!(input, sa)
    end
    return sa
end

function load_sa!(input::IO, sa::Vector{T}) where {T}
    # load a suffix array from the `input` into `sa`
    buf = Vector{UInt}(undef, 5)
    i = 0
    while !eof(input)
        read!(input, buf)
        value = T(0)
        @inbounds for j in 1:5
            value |= convert(T, buf[j]) << (8 * (j - 1))
        end
        sa[i+=1] = value
    end
    @assert i == length(sa)
    return sa
end


# other utils

function index_type(n)
    n -= 1
    n ≤ typemax(UInt8) ? UInt8 :
    n ≤ typemax(UInt16) ? UInt16 :
    n ≤ typemax(UInt32) ? UInt32 : UInt64
end

# suffix array sampling
function sample_sa(sa::Vector{T}, r) where {T}
    n = length(sa)
    samples = Vector{T}(undef, cld(n, r))
    sampled = falses(n)
    i′ = 0
    for i in 1:n
        @assert 0 ≤ sa[i] ≤ n - 1
        if sa[i] % r == 0
            samples[i′+=1] = sa[i]
            sampled[i] = true
        end
    end
    return samples, sampled
end

# Burrows-Wheeler Transform
function make_bwt(seq, sa)
    n = length(seq)
    @assert length(sa) == n
    ret = Vector{UInt}(undef, n)
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

struct LocationIterator{w,T}
    range::UnitRange{Int}
    index::FMIndex{w,T}
end

Base.length(iter::LocationIterator) = length(iter.range)

function Base.iterate(iter::LocationIterator, i::Int=1)
    if i > length(iter)
        return nothing
    end
    return sa_value(iter.range[i], iter.index) + 1, i + 1
end

end # FMIndexes
