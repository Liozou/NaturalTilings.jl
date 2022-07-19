
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
    ((_, ofs), x2), fst = findmin(Base.Fix1(getindex, kp), ering)
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
function identify_junction(r1::Vector{PeriodicVertex{D}}, r2, r3=nothing) where D
    len1 = length(r1)
    len2 = length(r2)
    start1 = 1
    _start2 = findfirst(==(r1[1]), r2)
    local start2::Int
    if _start2 isa Int
        start2 = _start2
    else
        start1 = 1 + length(r1)÷2
        _start2 = findfirst(==(r1[start1]), r2)
        if r3 isa Nothing && _start2 isa Nothing
            return nothing # more than 2 junctions
        else
            start2 = _start2
        end
    end
    after1 = r1[start1+1]
    after2 = r2[mod1(start2+1, len2)]
    direct = 2*(after1 == after2 || r1[mod1(start1-1, len1)] == r2[mod1(start2-1, len2)])-1
    _ia1 = mod1(start1 + 1, len1)
    _ia2 = mod1(start2 + direct, len2)
    _ib1 = mod1(start1 - 1, len1)
    _ib2 = mod1(start2 - direct, len2)
    while r1[_ia1] == r2[_ia2]
        _ia1 = mod1(_ia1 + 1, len1)
        _ia2 = mod1(_ia2 + direct, len2)
    end
    while r1[_ib1] == r2[_ib2]
        _ib1 = mod1(_ib1 - 1, len1)
        _ib2 = mod1(_ib2 - direct, len2)
    end

    ia1 = mod1(_ia1 - 1, len1)
    ib1 = mod1(_ib1 + 1, len1)
    ia2 = mod1(_ia2 - direct, len2)
    ib2 = mod1(_ib2 + direct, len2)
    iabs = [minmax(ia1, ib1)..., minmax(ia2, ib2)...]
    (xa, ofsa), (xb, ofsb) = minmax(r1[ia1], r1[ib1])
    @assert (PeriodicVertex(xa, ofsa), PeriodicVertex(xb, ofsb)) == minmax(r2[ia2], r2[ib2])
    if r3 isa Nothing
        set2 = Set{PeriodicVertex{D}}()
        while _ia2 != _ib2
            push!(set2, r2[_ia2])
            _ia2 = mod1(_ia2 + direct, len2)
        end
        while _ia1 != _ib1
            r1[_ia1] ∈ set2 && return nothing # more than 2 junctions
            _ia1 = mod1(_ia1 + 1, len1)
        end
    else
        ia3 = findfirst(==(r1[ia1]), r3)
        δ = ib1 < ia1 ? len1 - ia1 + ib1 : ib1 - ia1
        ib3 = mod1(ia3 + δ, length(r3))
        if r3[ib3] != r1[ib1]
            ib3 = mod1(ia3 - δ, length(r3))
            @assert r3[ib3] == r1[ib1]
        end
        ia3, ib3 = minmax(ia3, ib3)
        push!(iabs, ia3, ib3)
    end
    return PeriodicEdge{D}(xa, xb, ofsb .- ofsa), iabs
end

function delete_trivial_junctions!(junctions, jbitsets, junctions_per_ring, rings_with_junctions)
    delete_junctions = Int[]
    junctions_map = zeros(length(junctions))
    count = 1
    for i in 1:length(junctions)
        if length(jbitsets[i]) ≤ 3
            push!(delete_junctions, i)
        else
            junctions_map[i] = count
            count += 1
        end
    end
    deleteat!(junctions, delete_junctions)

    delete_rings = Int[]
    for (idx_k, k) in enumerate(rings_with_junctions)
        empty!(delete_junctions)
        js = junctions_per_ring[k]
        for (j, (x, i1, i2)) in enumerate(js)
            m = junctions_map[x]
            if m == 0
                push!(delete_junctions, j)
            else
                js[j] = (m, i1, i2)
            end
        end
        if length(delete_junctions) == length(js)
            push!(delete_rings, idx_k)
        else
            deleteat!(js, delete_junctions)
        end
    end
    deleteat!(rings_with_junctions, delete_rings)

    return junctions, junctions_per_ring, rings_with_junctions
end

function clean_junctions!(junctions::Vector{PeriodicEdge{D}}, rings, m, n) where D
    junctiongraph = [Tuple{Int, PeriodicVertex{D}}[] for _ in 1:n]
    for (idx_junction, (j_src, j_dst)) in enumerate(junctions)
        push!(junctiongraph[j_src], (idx_junction, j_dst))
    end
    ring_per_junctions = BitSet[BitSet() for _ in eachindex(junctions)]
    junctions_per_ring = Vector{Vector{Tuple{Int,Int,Int}}}(undef, m)
    rings_with_junctions = Int[]
    for k in 1:m
        ring = rings[k]
        vs = Dict{PeriodicVertex{D},Int}(w => j for (j, w) in enumerate(ring))
        flag = false
        for (i, (v, v_ofs)) in enumerate(ring)
            for (idx_junction, (_x, _x_ofs)) in junctiongraph[v]
                x = PeriodicVertex(_x, v_ofs + _x_ofs)
                j = get(vs, x, 0)
                j == 0 && continue
                flag = true
                if !isassigned(junctions_per_ring, k)
                    junctions_per_ring[k] = Tuple{Int,Int,Int}[(idx_junction, minmax(i, j)...)]
                else
                    push!(junctions_per_ring[k], (idx_junction, minmax(i, j)...))
                end
                push!(ring_per_junctions[idx_junction], k)
            end
        end
        flag && push!(rings_with_junctions, k)
    end

    return delete_trivial_junctions!(junctions, ring_per_junctions, junctions_per_ring, rings_with_junctions)
end

"""
    find_phantomedges(erings::Vector{Vector{Int}}, rings::Vector{Vector{PeriodicVertex{D}}}, kp::EdgeDict{D}, n, rings_idx) where D

Identify the set of junction points that split conjoined rings, i.e. two rings whose sum is
another cycle which is smaller or equal in size. One of the two rings must be have its
index in `rings_idx` to be considered. `n` is the number of vertices of the graph.

Return a triplet `(junctions, junctions_per_ring, keep)` where
- `junctions` is a list of `PeriodicEdge{D}` obtained by [`NaturalTilings.identify_junction`](@ref).
- `junctions_per_ring[i]` where ``i ∈ keep`` is a list of triplets `(x, j1, j2)` where `x`
  is the index of the corresponding junction in `junctions`, and `rings[i][j1]` and
  `rings[i][j2]` are the corresponding two junction points in the ring.
"""
function find_phantomedges(erings::Vector{Vector{Int}}, rings::Vector{Vector{PeriodicVertex{D}}}, kp::EdgeDict{D}, n, rings_idx) where D
    m = length(erings)
    eringdict = Dict{Vector{Int},Int}(x => i for (i,x) in enumerate(erings))
    junctions = PeriodicEdge{D}[] # the list of newly added edges
    junctions_dict = Dict{PeriodicEdge{D},Int}() # reverse map to junctions
    known_erings = Dict{Vector{Int},Int}()
    for (i,e) in enumerate(erings)
        e2 = copy(e)
        _start_ofs = canonical_ering!(e2, kp)
        @assert iszero(_start_ofs) # otherwise, the offset computed below at the next canonical_ering! call is wrong
        known_erings[e2] = i
    end
    buffer = Int[]
    eras = ERingAttributions(erings, kp)
    encountered = BitVector(undef, m)
    for i1 in rings_idx
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
                i3 = get(known_erings, buffer, nothing)
                _, i2, _, ofs_ref = find_ofs_ref(eri, _i2)
                r2 = PeriodicGraphs.OffsetVertexIterator{D}(ofs_ref, rings[i2])
                is = [i1, i2]
                if i3 isa Nothing
                    possible_junction = identify_junction(r1, r2)
                    possible_junction isa Nothing && continue
                    junction, iabs = possible_junction
                else
                    r3 = PeriodicGraphs.OffsetVertexIterator{D}(ofs_ear, rings[i3])
                    junction, iabs = identify_junction(r1, r2, r3)
                    push!(is, i3)
                end
                junction_idx = get!(junctions_dict, junction, length(junctions)+1)
                junction_idx > length(junctions) && push!(junctions, junction)
            end
        end
    end

    return clean_junctions!(junctions, rings, m, n)
end

function add_phantomedges!(erings::Vector{Vector{Int}}, rings::Vector{Vector{PeriodicVertex{D}}}, kp::EdgeDict{D}, new_rings_idx::Vector{Int}, junctions, junctions_per_ring, rings_with_junctions) where D
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
        for (i, src, dst) in junctions_per_ring[k]
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
    unsafe_copyto!(new_rings_idx, 1, invperm(I), length(I) - n + 1, n)
    permute!(rings, I)
    permute!(erings, I)
    return max_realedge
end


function explore_around_phantomedge!(tile::Vector{PeriodicVertex{D}}, start, tiling::Tiling{D}, max_realedge) where D
    tiledict = Dict(i => Set{SVector{D,Int}}() for (i,_) in tile)
    for (i, ofs) in tile
        push!(tiledict[i], ofs)
    end
    neighborindices = Vector{Tuple{Int,Int}}[]
    Q = [tile[start]]
    visited_rings = Set(Q)
    marked_edges = Dict{Int,Int}()
    for (i, (idx_ring, ofs_ring)) in enumerate(Q)
        ring = tile.rings[idx_ring]
        last_x = last(ring)
        neighboridx = [(0,0) for _ in 1:length(ring)]
        for (j, x) in enumerate(ring)
            x1, x2 = minmax(last_x, x)
            e = kp[(x1, x2)]
            e > max_realedge || continue
            (xa, ofsa), (xb, ofsb) = x1, x2
            encountered = get(marked_edges, e, nothing)
            if encountered isa Tuple{Int,Int}
                (idx_ring_encountered, idx_in_encountered) = encountered
                neighboridx[j] = encountered
                neighborindices[idx_ring_encountered][idx_in_encountered] = (i, j)
                continue
            end
            edge = PeriodicEdge{D}(xa, xb, ofsb - ofsa)
            flag = false
            for (idx2, ofs2) in tiling.ringsofedge[edge]
                ofss = get(tiledict, idx2, nothing)
                ofss isa Nothing && continue
                idx2 == idx_ring && ofsa == ofs2 && continue
                new_ofs = ofs_ring + ofsa - ofs2
                newofs in ofss || continue
                marked_edges[e] = (i, j)
                new_tovisit = PeriodicVertex{D}(idx2, newofs)
                new_tovisit in visited_rings || push!(Q, new_tovisit)
                flag = true
                break # there can only be one such other ring
            end
            @assert flag
            last_x = x
        end
        push!(neighborindices, neighboridx)
    end
    return neighborindices, Q
end

function remove_phantomedges(tiling::Tiling{D,T}, max_realedge::Int) where {D,T}
    # removerings_list = Int[]
    # removerings = zeros(Int, length(tiling.rings))
    # actualedges = Tuple{PeriodicEdge{D},SVector{D,Int}}[]
    # for (i, ering) in enumerate(tiling.erings)
    #     for e in ering
    #         removerings[e] == 0 || continue
    #         if e > max_realedge
    #             push!(removerings, i)
    #             removerings[e] = length(removerings_list)
    #             (x1, ofs1), (x2, ofs2) = kp[e]
    #             push!(actualedges, (PeriodicEdge{D}(x1, x2, ofs2 .- ofs1), ofs1))
    #             break
    #         end
    #     end
    # end
    # removerings[removerings_list] .= 1:length(removerings_list)
    for (i, tile) in enumerate(tiling.tiles)
        marked_edges .= 0
        for start in 1:length(tile)
            # for each ring of the cycle, find the largest set of rings of the tile
            # connected to it via phantom edges.
            neighborindices, Q = explore_around_phantomedge!(tile, start, tiling, max_realedge)

        end
    end
    
end
