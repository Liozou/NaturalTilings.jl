struct _ERingIncluding{D,T} <: AbstractVector{Vector{Int}}
    eras::T
    i::Int
    buffer::Vector{Int}
end

function ref_vertexpair(kp::EdgeDict{D}, x::Int)
    (v1, ofs1), (v2, ofs2) = kp[x]
    iszero(ofs1) && return x, ofs1
    return get!(kp, (PeriodicVertex{D}(v1), PeriodicVertex{D}(v2, ofs2 .- ofs1))), ofs1
end

struct ERingAttributions{D} <: AbstractVector{_ERingIncluding{D,ERingAttributions{D}}}
    erings::Vector{Vector{Int}}
    attrs::Vector{Vector{Tuple{Int,Int}}}
    kp::EdgeDict{D}

    function ERingAttributions(erings, kp::EdgeDict{D}) where D
        n = length(kp.direct)
        attrs = [Tuple{Int,Int}[] for _ in 1:n]
        for (i, er) in enumerate(erings)
            for (j, x) in enumerate(er)
                ref, = ref_vertexpair(kp, x)
                @assert ref ≤ n
                push!(attrs[ref], (i, j))
            end
        end
        return new{D}(erings, attrs, kp)
    end
end

Base.@propagate_inbounds function Base.getindex(eras::ERingAttributions, i::Integer)
    @boundscheck checkbounds(eras.attrs, i)
    ERingIncluding(eras, i)
end
Base.size(eras::ERingAttributions) = size(eras.attrs)
Base.IndexStyle(::Type{ERingAttributions{D}}) where {D} = Base.IndexLinear()

const ERingIncluding{D} = _ERingIncluding{D,ERingAttributions{D}}
ERingIncluding(eras::ERingAttributions{D}, i) where {D} = ERingIncluding{D}(eras, i)

function Base.getindex(eri::ERingIncluding{D}, j::Integer) where {D}
    newring_idx, idx = eri.ras.attrs[eri.i][j]
    newring = ri.ras.rings[newring_idx]
    ofs = newring[idx].ofs
    return OffsetVertexIterator{D}(.-ofs, newring)
end
Base.size(eri::ERingIncluding) = size(eri.ras.attrs[eri.i])
Base.IndexStyle(::Type{ERingIncluding{D}}) where {D} = Base.IndexLinear()



function add_phantomedges!(rings::Vector{Vector{PeriodicVertex{D}}}, erings, kp) where D
    n = length(rings)
    eringdict = Dict{Vector{Int},Int}(x => i for (i,x) in erings)
    junctions = PeriodicEdge{D}[] # the list of newly added edges
    junctions_per_ring = [Tuple{PeriodicVertex{D},Int,Int}[] for _ in 1:n]
    # junctions_per_ring[i] is the list of tuples (x, j, k) indicating that junction x is
    # present between vertices j and k of rings[i].
    era = ERingAttributions(erings, kp)
    for (i1, er1) in enumerate(erings)
        
        for i2 in i1:length(erings)
            er2 = erings[i2]
            er3 = sum_
        end
    end
end
