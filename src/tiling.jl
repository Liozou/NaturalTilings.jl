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
            push!(doubleedges, get!(tiling.kp, minmax(lst, x)))
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

function neighboring_edges(_edges::Vector{Int}, vertices::Vector{PeriodicVertex{D}}, tiling::Tiling{D,T}) where {D,T}
    edges = Set{Int}(_edges)
    returned_edges = Int[]
    for x in vertices
        for y in neighbors(tiling.g, x)
            e = get!(tiling.kp, minmax(x, y))
            e in edges && continue
            push!(returned_edges, e)
        end
    end
    sort!(returned_edges)
    unique!(returned_edges)
    return returned_edges
end

function lineof(e::Int, tiling)
    x1, x2 = tiling.kp[e]
    return shortline(tiling.pge[x1], tiling.pge[x2])
end

"""
    canonical_tile!(tiling::Tiling{D}, tile::Vector{PeriodicVertex{D}}, checkgeometry) where D

From a `tile` given as a list of its rings, return the position of the tile in `tiling`
and update `tiling` to add the new tile if it is absent.

The input `tile` may also be modified by this function.
"""
function canonical_tile!(tiling::Tiling{D}, tile::Vector{PeriodicVertex{D}}, checkgeometry) where D
    normalized, ofs = normalize_cycle!(tile)
    uniqueedges = unique_edges(normalized, tiling)
    isempty(uniqueedges) && return PeriodicVertex{D}(0)
    n = length(tiling.tiles) + 1
    idx = get!(tiling.tiledict, uniqueedges, n)
    if idx == n
        uniquevertices = unique_vertices(normalized, tiling)
        if checkgeometry
            for e in neighboring_edges(uniqueedges, uniquevertices, tiling)
                line = lineof(e, tiling)
                # tiling.kp[e] == (PeriodicVertex3D(1), PeriodicVertex3D(23)) && length(tile) == 20 && @show tile
                # if tile == PeriodicVertex3D[(3, (0,0,0)), (4, (0,0,0)), (6, (1,0,0)), (8, (1,0,0)), (9, (0,1,0)), (10, (0,1,0)), (11, (1,1,0)), (12, (1,1,0)), (13, (1,0,0)), (14, (1,0,0)), (23, (0,0,0)), (24, (0,0,0)), (45, (1,0,0)), (46, (1,1,0)), (48, (1,0,0)), (49, (0,0,0)), (50, (1,1,0)), (52, (0,1,0)), (54, (0,0,0)), (54, (0,1,0))] && e == 8316
                #     @show tiling.kp[e]
                #     @show intersect_polyedra(tiling, line, normalized)
                # end
                if intersect_polyedra(tiling, line, normalized) ≥ 2 # invalid tile
                    print(intersect_polyedra(tiling, line, normalized))
                    delete!(tiling.tiledict, uniqueedges)
                    return PeriodicVertex{D}(0)
                end
            end
        end
        checkgeometry && print("!")
        push!(tiling.tiles, normalized)
        push!(tiling.tileedges, uniqueedges)
        push!(tiling.tilevertices, uniquevertices)
        for (i,x) in enumerate(normalized)
            fst, lst = tiling.tilesofring[x.v]
            new = PeriodicVertex{D}(idx, x.ofs + ofs)
            if iszero(first(fst))
                tiling.tilesofring[x.v] = (new, lst)
            elseif iszero(first(lst))
                tiling.tilesofring[x.v] = (fst, new)
            else
                # @error "Ring $x of new tile $n joins more than 2 tiles."

                @show x, fst, lst

                # pop!(tiling.tiles)
                # pop!(tiling.tileedges)
                # tiling.tiledict[uniqueedges] = -1
                # for j in 1:(i-1)
                #     y = normalized[j]
                #     fst, lst = tiling.tilesofring[x.v]
                #     if iszero(first(lst))
                #         tiling.tilesofring[x.v] = (lst, lst)
                #     else
                #         tiling.tilesofring[x.v] = (fst, PeriodicVertex{D}(0))
                #     end
                # end

                return PeriodicVertex{D}(-1) # flags an error
            end
        end
    elseif idx > 0 && tiling.tiles[idx] != normalized
        # two tiles share the same edges but are made from different rings
        @show tiling.tiles[idx]
        @show normalized
        @error "Colliding tiles"
        delete!(tiling.tiledict, uniqueedges)
        return PeriodicVertex{D}(-2)
    end
    return PeriodicVertex{D}(idx, ofs)
end

function add_rtile!(rt::Vector{PeriodicVertex{D}}, gauss, known_htiles, tiling, m, checkgeometry, ofs=nothing) where D
    normalize_cycle!(rt)
    t = sort!([hash_position(PeriodicVertex{D}(v, isnothing(ofs) ? o : o .+ ofs), m) for (v,o) in rt])
    t ∈ known_htiles && return 0
    push!(known_htiles, t)
    gaussian_elimination!(gauss, t) && return 0
    isnothing(ofs) || return 0
    tile, ofs = canonical_tile!(tiling, rt, checkgeometry)
    @assert iszero(ofs)
    return tile
end

function tilingof(pge::PeriodicGraphEmbedding{D}, depth::Integer=10, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(pge.g), dist::DistanceRecord=DistanceRecord(pge.g,depth), stophere=missing) where D
    g = copy(pge.g)
    _rings, symms, erings, kp = strong_erings(g, depth, symmetries, dist)
    rings = Vector{PeriodicVertex{D}}[[reverse_hash_position(x, g) for x in r] for r in _rings]

    max_realedges = Int[]
    # export_vtf("/tmp/sos_1.vtf", PeriodicGraphEmbedding(g, pge.pos, pge.cell), nothing, 5) # TODO: remove

    new_rings_idx = collect(1:length(rings))
    junctions, junctions_per_ring, keep = find_phantomedges(erings, rings, kp, nv(g), copy(new_rings_idx))
    for e in junctions; add_edge!(g, e) end
    # export_vtf("/tmp/sos_2.vtf", PeriodicGraphEmbedding(g, pge.pos, pge.cell), nothing, 5) # TODO: remove
    max_realedge = add_phantomedges!(erings, rings, kp, new_rings_idx, junctions, junctions_per_ring, keep)
    while !iszero(max_realedge)
        junctions, junctions_per_ring, keep = find_phantomedges(erings, rings, kp, nv(g), copy(new_rings_idx))
        for e in junctions; add_edge!(g, e) end
        # export_vtf("/tmp/sos_$max_realedge.vtf", PeriodicGraphEmbedding(g, pge.pos, pge.cell), nothing, 5) # TODO: remove
        max_realedge = add_phantomedges!(erings, rings, kp, new_rings_idx, junctions, junctions_per_ring, keep)
        push!(max_realedges, max_realedge)
    end

    if !isempty(max_realedges)
        symms = NoSymmetryGroup(length(rings)) # TODO: adapt the symmetry instead
    end
    tiling = Tiling(pge, g, rings, erings, kp)
    uniquesymms = unique(symms)
    m = length(rings)
    exploration = [TilingAroundCycle{D}(u) for u in uniquesymms]
    missing2tiles = collect(1:length(exploration))
    tilecounter = zeros(Int, m)
    etiles_gauss = IterativeGaussianEliminationLength()
    num_tiles = 0
    num_iter = 0
    known_htiles = Set{Vector{Int}}()
    checkgeometry = false
    while !isempty(missing2tiles) && num_iter < 3
        num_iter += 1
        new_rtiles = Vector{PeriodicVertex{D}}[] # each sublist is the list of rings of a tile
        for r in missing2tiles
            tac = exploration[r]
            _new_rtiles = explore_around_cycle!(tac, tiling)
            append!(new_rtiles, _new_rtiles)
        end
        isempty(new_rtiles) && continue
        sort!(new_rtiles; by=length)
        new_idx = length(tiling.tiles) + 1
        for rt in new_rtiles
            tile_idx = add_rtile!(rt, etiles_gauss, known_htiles, tiling, m, stophere isa Bool ? (stophere || checkgeometry) : checkgeometry)
            if tile_idx < 0 # error occured: one ring joins 2 or more tiles
                stophere isa Bool && (stophere = checkgeometry = true; break)
                empty!(new_rtiles)
                break
            end
            if tile_idx == num_tiles + 1 # new tile
                num_tiles = tile_idx
                for symm in symms
                    symmrt = [symm()] # FIXME: this part is invalid, try with any symmetry
                    tile_idx2 = add_rtile!([symm(x) for x in rt], etiles_gauss, known_htiles, tiling, m, stophere isa Bool ? (stophere || checkgeometry) : checkgeometry)
                    num_tiles = max(num_tiles, tile_idx2)
                end
                for tile in @view tiling.tiles[tile_idx:end]
                    for outer_ofs in PeriodicGraphs.cages_around(g, 2)
                        add_rtile!(tile, etiles_gauss, known_htiles, tiling, m, checkgeometry, outer_ofs)
                    end
                end
            end
        end
        stophere isa Bool && stophere && checkgeometry && break

        @assert num_tiles == length(tiling.tiles) || isempty(new_rtiles)
        if !isempty(new_rtiles)
            for i in new_idx:num_tiles
                for (v, _) in tiling.tiles[i]
                    if tilecounter[v] == 2 # # error occured: one ring joins 2 or more tiles
                        stophere isa Bool && (stophere = checkgeometry = true; break)
                        empty!(new_rtiles)
                        break
                    end
                    tilecounter[v] += 1
                end
                stophere isa Bool && stophere && checkgeometry && break
                isempty(new_rtiles) && break
            end
        end
        stophere isa Bool && stophere && checkgeometry && break

        if isempty(new_rtiles) # error occured: one ring joins 2 or more tiles
            empty!(tiling.tiles)
            empty!(tiling.tileedges)
            empty!(tiling.tilevertices)
            empty!(tiling.tiledict)
            for i in 1:m
                tiling.tilesofring[i] = (PeriodicVertex{D}(0),PeriodicVertex{D}(0))
            end
            (D != 3 || checkgeometry) && return tiling # abandon and return empty tiling
            # otherwise, restart everything but now checking the geometry of each tile
            checkgeometry = true
            exploration = [TilingAroundCycle{D}(u) for u in uniquesymms]
            missing2tiles = collect(1:length(exploration))
            tilecounter = zeros(Int, m)
            etiles_gauss = IterativeGaussianEliminationLength()
            num_tiles = 0
            num_iter = 0
            known_htiles = Set{Vector{Int}}()
            continue
        end

        missing2tiles = [i for i in uniquesymms if tilecounter[i] < 2]
    end
    isempty(missing2tiles) || isempty(tiling.tiles) || @error "Possibly missed some large tiles"

    return tiling
end

tilingof(pge::PeriodicGraphEmbedding{D}, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(pge.g,10)) where {D} = tilingof(pge, 10, symmetries, dist)

# for debugging purpose # TODO: remove
tilingof(pge::PeriodicGraphEmbedding, debug::Bool) = tilingof(pge, 10, NoSymmetryGroup(pge.g), DistanceRecord(pge.g,10), debug)
