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
        uniquevertices = unique_vertices(normalized, tiling)
        push!(tiling.tilevertices, uniquevertices)
        for (i,x) in enumerate(normalized)
            fst, lst = tiling.tilesofring[x.v]
            new = PeriodicVertex{D}(idx, x.ofs + ofs)
            if iszero(first(fst))
                tiling.tilesofring[x.v] = (new, fst)
            elseif iszero(first(lst))
                tiling.tilesofring[x.v] = (fst, new)
            else
                # @error "Ring $x of new tile $n joins more than 2 tiles."

                # @show fst
                # @show lst
                # @show new
                # @show x
                # @show normalized
                pop!(tiling.tiles)
                pop!(tiling.tileedges)
                tiling.tiledict[uniqueedges] = -1
                for j in 1:(i-1)
                    y = normalized[j]
                    fst, lst = tiling.tilesofring[x.v]
                    if iszero(first(lst))
                        tiling.tilesofring[x.v] = (lst, lst)
                    else
                        tiling.tilesofring[x.v] = (fst, PeriodicVertex{D}(0))
                    end
                end
                return PeriodicVertex{D}(-1)
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
    max_realedge = add_phantomedges!(erings, rings, kp)
    if max_realedge != 0
        symms = NoSymmetryGroup(length(rings)) # TODO: adapt the symmetry instead
    end
    tiling = Tiling{D}(rings, erings, kp)
    uniquesymms = unique(symms)
    m = length(rings)
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
            # tile_idx < 0 && (println(-tile_idx); continue)
            # tile_idx < 0 && (println(-tile_idx, ' ', name); return nothing)
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
                    # println("X ", name)
                    # return nothing
                    @error "Ring $v joins more than 2 tiles."
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
