
"""
    ref_vertexpair(kp::EdgeDict{D}, x::Int) where D

Return `(y, ofs)` where `y::Int` refers to the same edge as `x` offset by `-ofs`, such that
all representatives of `x` yield the same `y`.
"""
function ref_vertexpair(kp::EdgeDict{D}, x::Int) where D
    (v1, ofs1), (v2, ofs2) = kp[x]
    iszero(ofs1) && return x, ofs1
    return kp[(PeriodicVertex{D}(v1), PeriodicVertex{D}(v2, ofs2 .- ofs1))], ofs1
end

struct ERingAttributions{D}
    erings::Vector{Vector{Int}}
    attrs::Vector{Vector{Tuple{Int,Int}}}
    kp::EdgeDict{D}

    function ERingAttributions(erings, kp::EdgeDict{D}) where D
        n = length(kp.direct)
        attributed = falses(n)
        lastattributed = 0
        attrs = Vector{Vector{Tuple{Int,Int}}}(undef, n)
        for (i, er) in enumerate(erings)
            for (j, x) in enumerate(er)
                ref, = ref_vertexpair(kp, x)
                @assert ref ≤ n
                if !attributed[ref]
                    attributed[ref] = true
                    attrs[ref] = [(i, j)]
                    lastattributed = max(lastattributed, ref)
                else
                    push!(attrs[ref], (i, j))
                end
            end
        end
        resize!(attrs, lastattributed)
        return new{D}(erings, attrs, kp)
    end
end

struct ERingIncludingExcluding{D} <: AbstractVector{Vector{Int}}
    eras::ERingAttributions{D}
    i::Int
    ofs::SVector{D,Int}
    encountered::BitVector
    min::Int
    buffer::Vector{Int}
    function ERingIncludingExcluding(eras::ERingAttributions{D}, i, ofs, encountered, min) where {D}
        new{D}(eras, i, ofs, encountered, min, Int[])
    end
end

function find_ofs_ref(eri::ERingIncludingExcluding{D}, i::Int) where {D}
    newring_idx, idx = eri.eras.attrs[eri.i][i]
    newering = eri.eras.erings[newring_idx]
    k = newering[idx]
    k == i && return newering, newring_idx, true, zero(SVector{D,Int})
    kp = eri.eras.kp
    (_, ofs_ref2), = kp[k]
    return newering, newring_idx, false, eri.ofs .- ofs_ref2
end

function Base.getindex(eri::ERingIncludingExcluding{D}, i::Int) where {D}
    newering, idx, skip, ofs_ref = find_ofs_ref(eri, i)
    if idx < eri.min # only inspect each pair of rings once by asserting their order
        empty!(eri.buffer)
        return eri.buffer
    end
    m = length(eri.eras.erings)
    h = hash_position(PeriodicVertex(idx, ofs_ref), m)
    h > length(eri.encountered) && append!(eri.encountered, false for _ in 1:(h-length(eri.encountered)))
    if eri.encountered[h]
        # if this particular cycle has already been encountered, do not yield it again to
        # avoid duplicating computations: instead, return an empty cycle.
        empty!(eri.buffer)
        return eri.buffer
    end
    eri.encountered[h] = true
    n = length(newering)
    resize!(eri.buffer, n)
    if skip
        unsafe_copyto!(eri.buffer, 1, newering, 1, n)
    else
        kp = eri.eras.kp
        for (j, e) in enumerate(newering)
            (x1, ofs1), (x2, ofs2) = kp[e]
            eri.buffer[j] = get!(kp, (PeriodicVertex{D}(x1, ofs1 .+ ofs_ref), PeriodicVertex{D}(x2, ofs2 .+ ofs_ref)))
        end
        sort!(eri.buffer)
    end
    return eri.buffer
end
Base.size(eri::ERingIncludingExcluding) = size(eri.eras.attrs[eri.i])
Base.IndexStyle(::Type{ERingIncludingExcluding{D}}) where {D} = Base.IndexLinear()

Base.@propagate_inbounds function erings_including_excluding(eras::ERingAttributions, i, encountered, min)
    j, ofs = ref_vertexpair(eras.kp, i)
    @boundscheck checkbounds(eras.attrs, j)
    ERingIncludingExcluding(eras, j, ofs, encountered, min)
end


function canonical_ering!(ering::Vector{Int}, kp::EdgeDict{D}) where D
    lenc = length(ering)
    ((_, ofs), x2), fst = findmin(j -> kp[j], ering)
    if !iszero(ofs)
        for i in 1:lenc
            (v1, ofs1), (v2, ofs2) = kp[ering[i]]
            ering[i] = kp[(PeriodicVertex{D}(v1, ofs1 .- ofs), PeriodicVertex{D}(v2, ofs2 .- ofs))]
        end
        sort!(ering)
    end
    ofs
end

"""
    identify_junction(r1::Vector{PeriodicVertex{D}}, r2, r3) where D

Given two rings `r1` and `r2` and their sum `r3`, a ring smaller or equal in size, find the
positions of the two junctions, that is the two vertices common to all three rings.
"""
function identify_junction(r1::Vector{PeriodicVertex{D}}, r2, r3) where D
    len1 = length(r1)
    len2 = length(r2)
    start1 = 1
    _start2 = findfirst(==(r1[1]), r2)
    local start2::Int
    if _start2 isa Int
        start2 = _start2
    else
        start1 = 1 + length(r1)÷2
        start2 = findfirst(==(r1[start1]), r2)
    end
    after1 = r1[start1+1]
    after2 = r2[mod1(start2+1, len2)]
    direct = 2*(after1 == after2 || r1[mod1(start1-1, len1)] == r2[mod1(start2-1, len2)])-1
    ia1 = mod1(start1 + 1, len1)
    ia2 = mod1(start2 + direct, len2)
    ib1 = mod1(start1 - 1, len1)
    ib2 = mod1(start2 - direct, len2)
    while r1[ia1] == r2[ia2]
        ia1 = mod1(ia1 + 1, len1)
        ia2 = mod1(ia2 + direct, len2)
    end
    while r1[ib1] == r2[ib2]
        ib1 = mod1(ib1 - 1, len1)
        ib2 = mod1(ib2 - direct, len2)
    end
    ia1, ib1 = minmax(mod1(ia1 - 1, len1), mod1(ib1 + 1, len1))
    ia2, ib2 = minmax(mod1(ia2 - direct, len2), mod1(ib2 + direct, len2))
    (xa, ofsa), (xb, ofsb) = minmax(r1[ia1], r1[ib1])
    @assert (PeriodicVertex(xa, ofsa), PeriodicVertex(xb, ofsb)) == minmax(r2[ia2], r2[ib2])
    ia3 = findfirst(==(r1[ia1]), r3)
    ib3 = mod1(ia3 + abs(ib1 - ia1), length(r3))
    if r3[ib3] != r1[ib1]
        ib3 = mod1(ia3 + len1 - abs(ib1 - ia1), length(r3))
        @assert r3[ib3] == r1[ib1]
    end
    ia3, ib3 = minmax(ia3, ib3)
    return PeriodicEdge{D}(xa, xb, ofsb .- ofsa), ia1, ib1, ia2, ib2, ia3, ib3
end

function find_phantomedges(erings::Vector{Vector{Int}}, rings::Vector{Vector{PeriodicVertex{D}}}, kp::EdgeDict{D}, new_rings_idx) where D
    n = length(erings)
    eringdict = Dict{Vector{Int},Int}(x => i for (i,x) in enumerate(erings))
    junctions = PeriodicEdge{D}[] # the list of newly added edges
    junctions_dict = Dict{PeriodicEdge{D},Int}() # reverse map to junctions
    junctions_per_ring = Vector{Vector{Tuple{Int,Int,Int}}}(undef, n)
    has_junctions = falses(n)
    # junctions_per_ring[i] is the list of tuples (x, j, k) indicating that junction x is
    # present between vertices j and k of rings[i].
    known_erings = Dict{Vector{Int},Int}()
    for (i,e) in enumerate(erings)
        e2 = copy(e)
        _start_ofs = canonical_ering!(e2, kp)
        @assert iszero(_start_ofs) # otherwise, the offset computed below at the next canonical_ering! call is wrong
        known_erings[e2] = i
    end
    buffer = Int[]
    eras = ERingAttributions(erings, kp)
    encountered = BitVector(undef, n)
    for i1 in new_rings_idx
        er1 = erings[i1]
        len1 = length(er1)
        r1 = rings[i1]
        encountered .= false
        for x1 in er1
            eri = erings_including_excluding(eras, x1, encountered, i1)
            for (_i2, er2) in enumerate(eri)
                len2 = length(er2)
                abs(len2 - len1) > 1 && continue # includes the case where er2 is empty
                PeriodicGraphs.symdiff_cycles!(buffer, er1, er2)
                length(buffer) ≤ min(len1, len2) || continue
                isempty(buffer) && continue
                ofs_ear = canonical_ering!(buffer, kp)
                i3 = get(known_erings, buffer, 0)
                i3 == 0 && continue
                r3 = PeriodicGraphs.OffsetVertexIterator{D}(ofs_ear, rings[i3])
                _, i2, _, ofs_ref = find_ofs_ref(eri, _i2)
                r2 = PeriodicGraphs.OffsetVertexIterator{D}(ofs_ref, rings[i2])
                junction, i1a, i1b, i2a, i2b, i3a, i3b = identify_junction(r1, r2, r3)
                junction_idx = get!(junctions_dict, junction, length(junctions)+1)
                junction_idx > length(junctions) && push!(junctions, junction)
                for i in (i1, i2, i3)
                    if !has_junctions[i]
                        junctions_per_ring[i] = Tuple{PeriodicVertex{D},Int,Int}[]
                        has_junctions[i] = true
                    end
                end
                push!(junctions_per_ring[i1], (junction_idx, i1a, i1b))
                push!(junctions_per_ring[i2], (junction_idx, i2a, i2b))
                push!(junctions_per_ring[i3], (junction_idx, i3a, i3b))
            end
        end
    end
    keep = Int[]
    for k in 1:n
        if has_junctions[k]
            x = junctions_per_ring[k]
            sort!(x)
            unique!(x)
            push!(keep, k)
        end
    end
    return junctions, junctions_per_ring, keep
end

function add_phantomedges!(erings::Vector{Vector{Int}}, rings::Vector{Vector{PeriodicVertex{D}}}, kp::EdgeDict{D}, new_rings_idx::Vector{Int}) where D
    junctions, junctions_per_ring, rings_with_junctions = find_phantomedges(erings, rings, kp, new_rings_idx)
    isempty(junctions) && return 0
    max_realedge = length(kp.direct) # beyond are indices of new edges, which are not real ones.
    for (src, (dst, ofs)) in junctions
        for outer_ofs in PeriodicGraphs.cages_around(PeriodicGraph{D}(), 2)
            # update kp with the indices of the new edges
            get!(kp, (PeriodicVertex{D}(src, outer_ofs), PeriodicVertex{D}(dst, ofs .+ outer_ofs)))
        end
    end

    new_rings = Vector{PeriodicVertex{D}}[]
    new_rings_set = Set{Vector{PeriodicVertex{D}}}()
    new_erings = Vector{Int}[]
    for k in rings_with_junctions
        ring = rings[k]
        n = length(ring)
        g = PeriodicGraph{0}(n, PeriodicEdge{0}[(i, i+1) for i in 1:(n-1)])
        add_edge!(g, PeriodicEdge{0}(1, n))
        for (_, src, dst) in junctions_per_ring[k]
            add_edge!(g, PeriodicEdge{0}(src, dst))
        end
        subrings, = strong_rings(g, 60) # 60 is a gross overestimation but it should be costless here since ndims(g) == 0
        for _r in subrings
            r = first(normalize_cycle!(ring[_r]))
            if r ∉ new_rings_set
                push!(new_rings_set, r)
                push!(new_rings, r)
                push!(new_erings, convert_to_ering(r, kp))
            end
        end
    end
    deleteat!(rings, rings_with_junctions)
    deleteat!(erings, rings_with_junctions)
    append!(rings, new_rings)
    append!(erings, new_erings)
    I = sortperm(rings; by=length)
    n = length(new_rings)
    resize!(new_rings_idx, n)
    unsafe_copyto!(new_rings_idx, 1, invperm(I), length(I) - n, n)
    permute!(rings, I)
    permute!(erings, I)
    return max_realedge
end
