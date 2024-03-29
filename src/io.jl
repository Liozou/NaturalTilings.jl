export export_tile

function export_tile(file, pge::PeriodicGraphEmbedding3D{T}, tiling::Tiling, rtile::AbstractVector{PeriodicVertex3D}, type=nothing) where T
    if splitext(file)[2] != ".vtf"
        file = file * ".vtf"
    end
    mkpath(splitdir(file)[1])
    rings::Vector{OffsetVertexIterator{3}} = [OffsetVertexIterator{3}(ofs, tiling.rings[x]) for (x, ofs) in rtile]
    vertices::Vector{PeriodicVertex3D} = unique!(sort!(reduce(vcat, rings; init=PeriodicVertex3D[])))
    centre = sum(pge.pos[x] .+ ofs for (x,ofs) in vertices) / length(vertices)

    open(file, write=true) do f
        println(f, """
        ###############################
        # written by NaturalTilings.jl
        ###############################
        """)

        corres = Dict{PeriodicVertex3D,Int}()
        for (i,x) in enumerate(vertices)
            corres[x] = i-1
            v, _ = x
            print(f, "atom $(i-1) name $v")
            isnothing(type) ? println(f) : println(f, " type $type")
        end
        println(f)

        encounterededges = Set{Tuple{Int,Int}}()
        for r in rings
            prev_x = corres[last(r)]
            for x in r
                next_x = corres[x]
                mm = minmax(prev_x, next_x)
                if mm ∉ encounterededges
                    push!(encounterededges, mm)
                    println(f, "bond ", prev_x, ':', next_x)
                end
                prev_x = next_x
            end
        end

        ((_a, _b, _c), (_α, _β, _γ)), mat = cell_parameters(pge.cell)
        println(f, "pbc $_a $_b $_c $_α $_β $_γ\n")

        println(f, "ordered")
        for (v, ofs) in vertices
            _coord = ((T <: Rational ? widen.(pge.pos[v]) : pge.pos[v]) .+ ofs)
            coord = Float64.(mat * (0.8 .* _coord .+ 0.2 .* centre))
            join(f, round.(coord; digits=15), ' ')
            println(f)
        end
    end
    nothing
end

function export_tile(file, pge::PeriodicGraphEmbedding3D, tiling::Tiling{D}, (v,ofs)::PeriodicVertex{D}) where D
    rtile = OffsetVertexIterator{D}(tiling.tiles[v], ofs)
    export_tile(file, pge, tiling, rtile, v)
end
function export_tile(file, pge::PeriodicGraphEmbedding3D, tiling::Tiling, i::Integer)
    export_tile(file, pge, tiling, tiling.tiles[i], i)
end

function export_block(file, pge::PeriodicGraphEmbedding3D, tiling::Tiling, tiles::AbstractVector{<:Integer})
    mat = Float64.(pge.cell.mat)
    open(file, "w") do f
        println(f, length(tiles))
        for i in tiles
            μ = mat * mean(tiling.ringcenters[j] + ofs for (j, ofs) in tiling.tiles[i])
            σ = mean(norm(mat * pge[x] - μ) for x in tiling.tilevertices[i])
            a, b, c = μ
            println(f, a, ' ', b, ' ', c, ' ', σ)
        end
    end
    nothing
end

function PeriodicGraphEmbeddings.export_cgd(file, tiling::Tiling{D}) where D
    open(file, "w") do f
        println(f, "TILING")
        println(f, "  NAME ", splitext(basename(file))[1])
        ((a, b, c), (α, β, γ)) = first(cell_parameters(tiling.pge.cell))
        println(f, "  CELL ", a, ' ', b, ' ', c, ' ', α, ' ', β, ' ', γ)
        println(f, "  FACES")
        for face in tiling.rings
            print(f, "    ", length(face))
            for x in face
                print(f, "  ")
                join(f, tiling.pge[x], ' ')
            end
            println(f)
        end
        println(f, "END")
    end
end
