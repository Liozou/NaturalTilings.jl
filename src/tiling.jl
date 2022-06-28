export tilingof

"""
    unique_edges(tile::Vector{PeriodicVertex{D}}, tiling::Tiling{D}) where D

Return the list of sorted edge hashes defining the border of a `tile`.

Return an empty list if one of the edges appears in three or more rings of the `tile`.
"""
function unique_edges(tile::Vector{PeriodicVertex{D}}, tiling::Tiling{D}) where D
    doubleedges = Int[]
    for t in tile
        ring = tiling.rings[t.v]
        lst = PeriodicVertex{D}(last(ring).v, last(ring).ofs .+ t.ofs)
        for r in ring
            x = PeriodicVertex{D}(r.v, r.ofs .+ t.ofs)
            push!(doubleedges, tiling.kp[minmax(lst, x)])
            lst = x
        end
    end
    sort!(doubleedges)
    ret = doubleedges[1:2:end]
    bef = first(ret)
    for i in 2:length(ret)
        y = ret[i]
        @assert y == doubleedges[2*i]
        bef ≥ y && return Int[] # error: one edge is included more than twice
        bef = y
    end
    return ret
end

"""
    canonical_tile!(tiling::Tiling{D}, tile::Vector{PeriodicVertex{D}}) where D

From a `tile` given as a list of its rings, return the position of the tile in `tiling`
and update `tiling` to add the new tile if it is absent.

The input `tile` may also be modified by this function.
"""
function canonical_tile!(tiling::Tiling{D}, tile::Vector{PeriodicVertex{D}}) where D
    normalized, ofs = normalize_cycle!(tile)
    uniqueedges = unique_edges(normalized, tiling)
    isempty(uniqueedges) && return PeriodicVertex{D}(0)
    n = length(tiling.tiles) + 1
    idx = get!(tiling.tiledict, uniqueedges, n)
    if idx == n
        push!(tiling.tiles, normalized)
        push!(tiling.tileedges, uniqueedges)
        push!(tiling.tilevertices, unique_vertices(normalized, tiling))
        for x in normalized
            fst, lst = tiling.tilesofring[x.v]
            new = PeriodicVertex{D}(idx, x.ofs + ofs)
            if iszero(fst.v)
                tiling.tilesofring[x.v] = (new, lst)
            elseif iszero(lst.v)
                tiling.tilesofring[x.v] = (fst, new)
            else
                if !(fst == new || lst == new)
                    # @show fst
                    # @show lst
                    # @show new
                    # @show x
                    # @show normalized
                end
            end
        end
    elseif tiling.tiles[idx] != normalized
        # two tiles share the same edges but are made from different rings
        @error "Colliding tiles"
    end
    return PeriodicVertex{D}(idx, ofs)
end

function _detect_ring_crossing!(intersections, pos, r1, i1, r2, i2, opp, idx1, idx2)
    bef1 = r1[mod1(i1-1, length(r1))]
    bef2 = r2[mod1(i2-1, length(r2))]
    bef1 == bef2 && return false
    aft1 = r1[mod1(i1+1, length(r1))]
    aft1 == bef2 && return false
    aft2 = r2[mod1(i2+1, length(r2))]
    (aft1 == aft2 || bef1 == aft2) && return false

    p_ref = pos[first(r1[i1])] .+ last(r1[i1])
    p_opp = pos[first(opp)] .+ last(opp) .- p_ref
    p_bef1 = pos[first(bef1)] .+ last(bef1) .- p_ref
    p_bef2 = pos[first(bef2)] .+ last(bef2) .- p_ref
    p_aft1 = pos[first(aft1)] .+ last(aft1) .- p_ref
    p_aft2 = pos[first(aft2)] .+ last(aft2) .- p_ref

    λ = norm(p_opp)
    proj_bef1 = p_bef1 - dot(p_bef1, p_opp)/λ*p_opp
    proj_bef2 = p_bef1 - dot(p_bef2, p_opp)/λ*p_opp
    proj_aft1 = p_bef1 - dot(p_aft1, p_opp)/λ*p_opp
    proj_aft2 = p_bef1 - dot(p_aft2, p_opp)/λ*p_opp

    if proj_bef1 == -proj_aft1
        if proj_bef2 == -proj_aft2
            push!(intersections, minmax(idx1, idx2))
            return true
        end
        proj_bef1, proj_bef2 = proj_bef2, proj_bef1
        proj_aft1, proj_aft2 = proj_aft2, proj_aft1
    end

    bef1_bef2 = dot(proj_bef1, proj_bef2)
    aft1_bef2 = dot(proj_aft1, proj_bef2)
    bef1_aft2 = dot(proj_bef1, proj_aft2)
    aft1_aft2 = dot(proj_aft1, proj_aft2)

    if (bef1_bef2 ≥ 0 && aft1_bef2 ≥ 0 && (bef1_aft2 < 0 || aft1_aft2 < 0)) ||
       (bef1_aft2 ≥ 0 && aft1_aft2 ≥ 0 && (bef1_bef2 < 0 || aft1_bef2 < 0))
        push!(intersections, minmax(idx1, idx2))
        @show r1, r2
        return true
    end
    false
end

function _find_opp(r::PeriodicGraphs.OffsetVertexIterator{D}, pos) where D
    lenr = length(r)
    opp = mod1(pos + fld(lenr, 2), lenr)
    oppb = iseven(lenr) ? 0 : mod1(pos + cld(lenr, 2), lenr)
    oppositeb = iszero(oppb) ? PeriodicVertex{D}(0) : r[oppb]
    vi = r[pos]
    if oppb != 0 && oppositeb < vi
        oppb = 0
        oppositeb = PeriodicVertex{D}(0)
    end
    opposite = r[opp]
    if opposite < vi
        iszero(oppb) && return 0, 0, PeriodicVertex{D}(0), PeriodicVertex{D}(0)
        opp = oppb
        opposite = r[oppb]
        oppb = 0
        oppositeb = PeriodicVertex{D}(0)
    end
    return opp, oppb, opposite, oppositeb
end

function detect_ring_crossing(pge::PeriodicGraphEmbedding{D}, rings::Vector{Vector{PeriodicVertex{D}}}) where D
    ras = RingAttributions{D}(length(pge), rings)
    intersections = Set{Tuple{Int,Int}}()
    for i in 1:length(pge)
        ri = ras[i]
        attri = ras.attrs[i]
        vi = PeriodicVertex{D}(i)
        stored_opps = [_find_opp(ring, last(x)) for (ring, x) in zip(ri, attri)]
        for (j1, (r1, (idx1, i1))) in enumerate(zip(ri, attri))
            opp1, oppb1, opposite1, oppositeb1 = stored_opps[j1]
            iszero(opp1) && continue
            for j2 in (j1+1):length(ri)
                opp2, oppb2, opposite2, oppositeb2 = stored_opps[j2]
                iszero(opp2) && continue
                r2 = ri[j2]
                abs(length(r2) - length(r1)) ≤ 1 || continue
                idx2, i2 = attri[j2]
                if opposite1 == opposite2
                    _detect_ring_crossing!(intersections, pge.pos, r1, i1, r2, i2, opposite1, idx1, idx2) && continue
                    _detect_ring_crossing!(intersections, pge.pos, r1, opp1, r2, opp2, vi, idx1, idx2) && continue
                end
                if !iszero(oppb1)
                    if oppositeb1 == opposite2
                        _detect_ring_crossing!(intersections, pge.pos, r1, i1, r2, i2, oppositeb1, idx1, idx2) && continue
                        _detect_ring_crossing!(intersections, pge.pos, r1, oppb1, r2, opp2, vi, idx1, idx2) && continue
                    end 
                    if !iszero(oppb2) && oppositeb1 == oppositeb2
                        _detect_ring_crossing!(intersections, pge.pos, r1, i1, r2, i2, oppositeb1, idx1, idx2) && continue
                        _detect_ring_crossing!(intersections, pge.pos, r1, oppb1, r2, oppb2, vi, idx1, idx2) && continue
                    end
                end
                if !iszero(oppb2) && opposite1 == oppositeb2
                    _detect_ring_crossing!(intersections, pge.pos, r1, i1, r2, i2, opposite1, idx1, idx2) && continue
                    _detect_ring_crossing!(intersections, pge.pos, r1, opp1, r2, oppb2, vi, idx1, idx2) && continue
                end
            end
        end
    end
    @show intersections

end

function add_rtile!(rt::Vector{PeriodicVertex{D}}, gauss, known_htiles, tiling, m, ofs=nothing) where D
    normalize_cycle!(rt)
    t = sort!([hash_position(PeriodicVertex{D}(v, isnothing(ofs) ? o : o .+ ofs), m) for (v,o) in rt])
    t ∈ known_htiles && return 0
    push!(known_htiles, t)
    gaussian_elimination!(gauss, t) && return 0
    isnothing(ofs) || return 0
    tile, ofs = canonical_tile!(tiling, rt)
    @assert iszero(ofs)
    return tile
end

function tilingof(g::PeriodicGraph{D}, depth::Integer=10, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    _rings, symms, erings, kp = strong_erings(g, depth, symmetries, dist)
    rings = Vector{PeriodicVertex{D}}[[reverse_hash_position(x, g) for x in r] for r in _rings]
    tiling = Tiling{D}(rings, erings, kp)
    # TODO: include detect_ring_crossing
    uniquesymms = unique(symms)
    m = length(_rings)
    exploration = [TilingAroundCycle{D}(u) for u in uniquesymms]
    missing2tiles = collect(1:length(exploration))
    tilecounter = zeros(Int, m)
    etiles_gauss = IterativeGaussianEliminationLength()
    num_tiles = 0
    num_iter = 0
    known_htiles = Set{Vector{Int}}()
    while !isempty(missing2tiles) && num_iter < 5
        num_iter += 1
        new_rtiles = Vector{PeriodicVertex{D}}[] # each sublist is the list of rings of a tile
        for e in missing2tiles
            tac = exploration[e]
            append!(new_rtiles, explore_around_cycle!(tac, tiling))
        end
        sort!(new_rtiles; by=length)
        new_idx = length(tiling.tiles) + 1
        for rt in new_rtiles
            tile_idx = add_rtile!(rt, etiles_gauss, known_htiles, tiling, m)
            if tile_idx == num_tiles + 1 # new tile
                num_tiles = tile_idx
                for symm in symms
                    symmrt = [symm()] # FIXME: this part is invalid, try with any symmetry
                    tile_idx2 = add_rtile!([symm(x) for x in rt], etiles_gauss, known_htiles, tiling, m)
                    num_tiles = max(num_tiles, tile_idx2)
                end
                for tile in @view tiling.tiles[tile_idx:end]
                    for outer_ofs in PeriodicGraphs.cages_around(g, 2)
                        add_rtile!(tile, etiles_gauss, known_htiles, tiling, m, outer_ofs)
                    end
                end
            end
        end
        @assert num_tiles == length(tiling.tiles)
        for i in new_idx:num_tiles
            for (v, _) in tiling.tiles[i]
                if tilecounter[v] == 2
                    @info "Ring $v joins more than 2 tiles."
                end
                tilecounter[v] += 1
            end
        end
        missing2tiles = [i for i in uniquesymms if tilecounter[i] < 2]
    end
    isempty(missing2tiles) || @error "Possibly missed some large tiles"

    return tiling
end

tilingof(g::PeriodicGraph{D}, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(g,10)) where {D} = tilingof(g, 10, symmetries, dist)
