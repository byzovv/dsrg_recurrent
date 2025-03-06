using Glob
using DelimitedFiles
using LinearAlgebra


# Function for obtaining a list of files
function rdir(dir::AbstractString, pat::Glob.FilenameMatch)
    result = String[]
    for (root, dirs, files) in walkdir(dir)
        append!(result, filter!(f -> occursin(pat, f), joinpath.(root, files)))
    end
    return result
end
rdir(dir::AbstractString, pat::AbstractString) = rdir(dir, Glob.FilenameMatch(pat))


filenames = rdir("output", "*.txt")


# Take all the matrices of adjacency in the output catalog
# and check that the corresponding digrans are strongly regular with the desired parameters
for filename in filenames
    bname = basename(filename)
    bname_parts = split(bname, ['_', '.'])

    v = parse(Int16, bname_parts[1])
    k = parse(Int16, bname_parts[2])
    t = parse(Int16, bname_parts[3])
    λ = parse(Int16, bname_parts[4])
    μ = parse(Int16, bname_parts[5])

    A = readdlm(filename, Int16)
    J = ones(Int16, v, v)
    II = Matrix{Int16}(I, v, v)

    println("Checking dsrg($(v), $(k), $(t), $(λ), $(μ)):")

    println(A * A == t * II + λ * A + μ * (J - II - A))
    println(A * J == k * J)
    println(J * A == k * J)
    println("----------------------------------------------------------")

end
