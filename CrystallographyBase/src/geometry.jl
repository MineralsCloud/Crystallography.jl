export vertices, edge, edges, faces

const FACES = Base.vect(
    (1, 2, 3, 4), (5, 6, 7, 8), (1, 2, 6, 5), (3, 4, 8, 7), (2, 3, 7, 6), (5, 8, 4, 1)
)

function vertices(lattice::Lattice)
    O⃗ = zeros(eltype(lattice), 3)
    A⃗, B⃗, C⃗ = basisvectors(lattice)
    A⃗B⃗, A⃗C⃗, B⃗C⃗, A⃗B⃗C⃗ = A⃗ + B⃗, A⃗ + C⃗, B⃗ + C⃗, A⃗ + B⃗ + C⃗
    return Base.vect(O⃗, A⃗, B⃗, C⃗, A⃗B⃗, A⃗C⃗, B⃗C⃗, A⃗B⃗C⃗)
end

edge(A⃗, B⃗) = hcat(([Aᵢ, Bᵢ] for (Aᵢ, Bᵢ) in zip(A⃗, B⃗))...)

function edges(lattice::Lattice)
    O⃗, A⃗, B⃗, C⃗, A⃗B⃗, A⃗C⃗, B⃗C⃗, A⃗B⃗C⃗ = vertices(lattice)
    return Base.vect(
        edge(O⃗, A⃗),
        edge(O⃗, C⃗),
        edge(A⃗, A⃗C⃗),
        edge(C⃗, A⃗C⃗),
        edge(O⃗, B⃗),
        edge(A⃗, A⃗B⃗),
        edge(C⃗, B⃗C⃗),
        edge(A⃗C⃗, A⃗B⃗C⃗),
        edge(B⃗, A⃗B⃗),
        edge(B⃗, B⃗C⃗),
        edge(A⃗B⃗, A⃗B⃗C⃗),
        edge(B⃗C⃗, A⃗B⃗C⃗),
    )
end

function faces(lattice::Lattice)
    verts = vertices(lattice)
    return map(face -> [verts[i] for i in face], FACES)
end
