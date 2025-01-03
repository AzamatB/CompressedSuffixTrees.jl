module WaveletMatrices

export WaveletMatrix, getindex, rank, select, freq

using IndexableBitVectors
using IntArrays

import Base: lastindex, length, size, sizeof
import IndexableBitVectors: rank

# build an internal structure of the wavelet matrix
function build(::Type{B}, data::AbstractVector{T}, w::Int, destructive::Bool) where {T<:Unsigned,B}
    ivec = IntVector{w}(data)
    if destructive
        # free memory
        empty!(data)
    end
    bits = Vector{B}(undef, w)
    _build!(B, ivec, bits)
    return tuple(bits...)
end

function _build!(::Type{B}, ivec, bits) where {B}
    w = length(bits)
    len = length(ivec)
    # allocate working space
    bv = BitVector(undef, len)
    ivec′ = similar(ivec)
    for d in 1:w
        # scan d-th bit
        nzeros = 0
        for i in 1:len
            if (ivec[i] >> (w - d)) & 1 == 1
                bv[i] = true  # right
            else
                bv[i] = false # left
                nzeros += 1
            end
        end
        # stably sort integers by bv
        r = nzeros
        l = 0
        for i in 1:len
            if bv[i]
                ivec′[r+=1] = ivec[i]
            else
                ivec′[l+=1] = ivec[i]
            end
        end
        # store the bit vector and go next
        bits[d] = bv
        copy!(ivec, ivec′)
    end
end

# starting points for the rank operation can be precomputed
function locate_sp(a, bits, nzeros)
    w = length(bits)
    sp = 0
    for d in 1:w
        bit = (a >> (w - d)) & 1 == 1
        sp = rank(bit, bits[d], sp)
        if bit
            sp += nzeros[d]
        end
    end
    return sp
end

"""
Indexable sequence.

`WaveletMatrix` is a data structure like the wavelet tree, which supports fast rank/select
queries for a sequence of alphabets. Elements (typed `T`) are encoded in `w` bits and therefore
a sequence of unsigned integers can be stored in a space-efficient manner if `w` is small.
The actual memory space is determined according to the underlying bit vector type.
The default bit vector type is `SucVector`, which requires 5/4 bits per bit, so the total
size will be about `w * 5/4 * length` bits.

See (Claude et al, 2012, doi:10.1007/978-3-642-34109-0_18) for more details.
"""
struct WaveletMatrix{w,T<:Unsigned,B<:AbstractBitVector} <: AbstractVector{T}
    bits::NTuple{w,B}
    nzeros::NTuple{w,Int}
    sps::Vector{Int}
    function WaveletMatrix{w,T,B}(bits) where {w,T<:Unsigned,B<:AbstractBitVector}
        @assert 1 ≤ w ≤ sizeof(T) * 8 ≤ 64
        @assert length(bits) == w
        nzeros = Int[]
        for bv in bits
            @assert length(bits[1]) == length(bv)
            push!(nzeros, rank0(bv, length(bv)))
        end
        if w ≤ 16
            # size of lookup table ≤ 512KiB (= sizeof(Int) * 2^16)
            alphabetsize = 2^w
            sps = Vector{Int}(undef, alphabetsize)
            for a in 0:alphabetsize-1
                sps[a+1] = locate_sp(T(a), bits, nzeros)
            end
        else
            sps = Int[]
        end
        return new(tuple(bits...), tuple(nzeros...), sps)
    end
    function WaveletMatrix{w}(data::AbstractVector{T}; destructive::Bool=false) where {w,T<:Unsigned}
        bits = build(default_bitvector, data, w, destructive)
        return WaveletMatrix{w,T,default_bitvector}(bits)
    end
    function WaveletMatrix(data::AbstractVector{T}, w::Integer=sizeof(T) * 8) where {T<:Unsigned}
        return WaveletMatrix{w}(data)
    end
end

const default_bitvector = SucVector

length(wm::WaveletMatrix) = length(wm.bits[1])
size(wm::WaveletMatrix) = (length(wm),)

function sizeof(wm::WaveletMatrix{w}) where {w}
    s = 0
    for d in 1:w
        s += sizeof(wm.bits[d])
    end
    s += sizeof(wm.nzeros)
    s += sizeof(wm.sps)
    return s
end

@inline function Base.getindex(wm::WaveletMatrix{w,T}, i::Int) where {w,T}
    if i < 0 || lastindex(wm) < i
        throw(BoundsError(i))
    end
    ret = T(0)
    @inbounds begin
        for d in 1:w-1
            bits = wm.bits[d]
            bit = bits[i]
            ret = ret << 1 | bit
            if bit
                i = wm.nzeros[d] + rank1(bits, i)
            else
                i = rank0(bits, i)
            end
        end
        bits = wm.bits[w]
        bit = bits[i]
        ret = ret << 1 | bit
    end
    return ret
end

@inline function Base.getindex(wm::WaveletMatrix{w,T,SucVector}, i::Int) where {w,T}
    if i < 0 || lastindex(wm) < i
        throw(BoundsError(i))
    end
    ret = T(0)
    @inbounds begin
        for d in 1:w-1
            bits = wm.bits[d]
            bit, rnk1 = IndexableBitVectors.accrank1(bits, i)
            ret = ret << 1 | bit
            if bit
                i = wm.nzeros[d] + rnk1
            else
                i = i - rnk1
            end
        end
        bits = wm.bits[w]
        bit = bits[i]
        ret = ret << 1 | bit
    end
    return ret
end

@inline Base.getindex(wm::WaveletMatrix{w,T}, i::Integer) where {w,T} = getindex(wm, convert(Int, i))

function rank(a::Unsigned, wm::WaveletMatrix{w}, i::Int) where {w}
    i = clamp(i, 0, lastindex(wm))
    if i == 0
        return 0
    end

    if !isempty(wm.sps)
        # use precomputed sp
        sp = wm.sps[a+1]
        ep = i
        @inbounds for d in 1:w
            bit = (a >> (w - d)) & 1 == 1
            ep = rank(bit, wm.bits[d], ep)
            if bit
                ep += wm.nzeros[d]
            end
        end
    else
        sp, ep = 0, i
        # scan from the most significant bit to the least significant bit
        @inbounds for d in 1:w
            bit = (a >> (w - d)) & 1 == 1
            sp = rank(bit, wm.bits[d], sp)
            ep = rank(bit, wm.bits[d], ep)
            if bit
                nz = wm.nzeros[d]
                sp += nz
                ep += nz
            end
        end
    end
    return ep - sp
end

rank(a::Unsigned, wm::WaveletMatrix, i::Integer) = rank(a, wm, convert(Int, i))

function select(a::Unsigned, wm::WaveletMatrix, j::Int)
    if j ≤ 0
        return 0
    end
    # binary search: j ∈ (rank(l), rank(u)]
    l = 0
    u = length(wm)
    rank_u = rank(a, wm, u)
    while j ≤ rank_u
        m = div(l + u, 2)
        if l == m
            return u
        end
        rank_m = rank(a, wm, m)
        if j ≤ rank_m
            u = m
            rank_u = rank_m
        else
            l = m
        end
    end
    return 0
end

select(a::Unsigned, wm::WaveletMatrix, j::Integer) = select(a, wm, convert(Int, j))

function freq(a::Unsigned, wm::WaveletMatrix, i::Integer, j::Integer)
    return j < i ? 0 : rank(a, wm, j) - rank(a, wm, i - 1)
end

end # WaveletMatrices
