using JuMP, KNITRO
using MathOptInterface
const MOI = MathOptInterface
using DelimitedFiles
using BlockArrays
using LinearAlgebra
using Kronecker


# The parameters of the first graph in the sequence
# Note: only such parameters must be indicated,
# for which there are appropriate adjacency matrices in the input catalog
first_v = 8
first_k = 4
t = 3
λ = 1


s = t - λ

# Function corresponding to alpha_k operations
function getPart(matr, k)
    nrows = size(matr, 1)
    res = matr[1:(nrows÷k), 1:end]
    return res
end

# Function of forming a exchange matrix Kn
function getKn(n)
    res = zeros(Int16, n, n)
    for i=1:n
        res[i, n-i+1] = 1
    end
    return res
end

# Function of forming a matrix Pn
function getPn(n)
    res = ones(Int16, 2^n, 1) ⊗ getKn(2^n) ⊗ ones(Int16, t, t*2^n)
    return res
end

# Read the first adjacency matrix from the file
A1 = readdlm(joinpath("input", "$(lpad(first_v, 2, '0'))_$(first_k)_$(t)_$(λ)_$(t).txt"), Int16)

model = Model(KNITRO.Optimizer)

# For some sets of parameters, we change the Knitro library settings
if ((t == 5) && (λ == 4)) || ((t == 2) && (λ == 1)) || ((t == 4) && (λ == 3))
    set_attribute(model, "mip_method", 2)
else
    set_attribute(model, "mip_method", 1)
end
if ((t == 4) && (λ == 3))
    set_attribute(model, "mip_heuristic_maxit", 200)
end

@variable(model, 0 <= B1[1:first_v, 1:t*4] <= 1, Bin, start=1)
@variable(model, 0 <= C1[1:t*4, 1:first_v] <= 1, Bin, start=1)

P1 = getPn(1)

# Conditions for sums in rows and columns of matrices B1 and C1
for i=1:first_v
    @constraint(model, sum(B1[i, j] for j in 1:t*4) == t*2)
end
for i=1:t*4
    @constraint(model, sum(C1[i, j] for j in 1:first_v) == first_k)
end
for j=1:t*4
    @constraint(model, sum(B1[i, j] for i in 1:first_v) == first_k)
end
for j=1:first_v
    @constraint(model, sum(C1[i, j] for i in 1:t*4) == t*2)
end


# Conditions of blockness of matrices B1 and C1
# For some sets of parameters, we clearly do not indicate these conditions,
# but the program still finds such B1 and C1 for which these conditions are true
if !(((t == 6) && (λ == 5)) || ((t == 7) && (λ == 5)))
    for i=1:first_v
        @NLconstraint(model, (sum(B1[i, j] for j in 1:(t*2))) * (sum(B1[i, j] for j in 1:(t*2)) - t*2) == 0)
        @NLconstraint(model, (sum(B1[i, j] for j in (t*2+1):(t*4))) * (sum(B1[i, j] for j in (t*2+1):(t*4)) - t*2) == 0)
    end
    for j=1:first_v
        @constraint(model, sum(C1[i, j] for i in 1:t*2) == t)
        @constraint(model, sum(C1[i, j] for i in (t*2+1):(t*4)) == t)
    end
end


# Equation B1 * C1 = t*J
for i=1:first_v
    for j=1:first_v
        @NLconstraint(model, sum(B1[i, k] * C1[k, j] for k in 1:t*4) == t)
    end
end

# Equation A1 * B1 + s*B1 = t*J
for i=1:first_v
    for j=1:t*4
        @constraint(model, sum(A1[i, k] * B1[k, j] for k in 1:first_v) + s * B1[i, j] == t)
    end
end

# Equation B1 * P1 = t*J
for i=1:first_v
    for j=1:t*4
        @constraint(model, sum(B1[i, k] * P1[k, j] for k in 1:t*4) == t)
    end
end

# Equation P1 * C1 = t*J
for i=1:t*4
    for j=1:first_v
        @constraint(model, sum(P1[i, k] * C1[k, j] for k in 1:t*4) == t)
    end
end

# Equation C1 * A1 + s*C1 = t*J
for i=1:t*4
    for j=1:first_v
        @constraint(model, sum(C1[i, k] * A1[k, j] for k in 1:first_v) + s * C1[i, j] == t)
    end
end

# Equation C1 * B1 + s*P1 = t*J
for i=1:t*4
    for j=1:t*4
        @NLconstraint(model, sum(C1[i, k] * B1[k, j] for k in 1:first_v) + s * P1[i, j] == t)
    end
end


# Optimize the formal function F = 1
@objective(model, Min, 1)
JuMP.optimize!(model)

# If you could not find matrices B1 and C1, we complete the work of the program
gap = MOI.get(model, MOI.RelativeGap())
if gap > 1e-6
    print("Unable to find matrices B1 and C1! The program terminates!")
    exit()
end

B1_values = zeros(Int16, first_v, t*4)
C1_values = zeros(Int16, t*4, first_v)

for i=1:first_v
    for j=1:t*4
        B1_values[i, j] = round(value(B1[i, j]))
    end
end

for i=1:t*4
    for j=1:first_v
        C1_values[i, j] = round(value(C1[i, j]))
    end
end

# Form the matrix A2 from blocks
A2 = mortar(reshape([A1, zero(A1), zero(C1_values), C1_values,
                     zero(A1), A1, C1_values, zero(C1_values),
                     B1_values, zero(B1_values), zero(P1), P1,
                     zero(B1_values), B1_values, P1, zero(P1)], 4, 4))

second_v = 2*first_v + 2*t*4
second_k = first_k + t*2

# Save A2 to the file
mkpath(joinpath("output", "$(t)_$(λ)_$(t)"))
writedlm(joinpath("output", "$(t)_$(λ)_$(t)", "$(second_v)_$(second_k)_$(t)_$(λ)_$(t).txt"), A2, " ")


# The function of forming a matrix Bn
function getBn(n)
    if n == 1
        res = B1_values
    else
        res = mortar(reshape([getKn(2) ⊗ getBn(n-1) ⊗ ones(Int16, 1, 2),
                              Matrix{Int16}(I, 2, 2) ⊗ getPn(n-1) ⊗ ones(Int16, 1, 2)], 2, 1))
    end
    return res
end

# The function of forming the matrix Cn
function getCn(n)
    if n == 1
        res = C1_values
    else
        res = mortar(reshape([ones(Int16, 2^n, 1) ⊗ Matrix{Int16}(I, 2, 2) ⊗ getPart(getCn(n-1), 2^(n-1)),
                              ones(Int16, 2^n, 1) ⊗ Matrix{Int16}(I, 2, 2) ⊗ getPart(getPn(n-1), 2^(n-1))], 1, 2))
    end
    return res
end

# The function of forming a matrix An
function getAn(n)
    if n == 1
        res = A1
    elseif n == 1
        res = A2
    else
        res = mortar(reshape([getAn(n-1), zero(getAn(n-1)), zero(getCn(n-1)), getCn(n-1),
                              zero(getAn(n-1)), getAn(n-1), getCn(n-1), zero(getCn(n-1)),
                              getBn(n-1), zero(getBn(n-1)), zero(getPn(n-1)), getPn(n-1),
                              zero(getBn(n-1)), getBn(n-1), getPn(n-1), zero(getPn(n-1))], 4, 4))
    end
    return res
end

# We form four more elements of the sequence: A3, A4, A5, A6
for n=3:6
    An = getAn(n)

    cur_v = size(An, 1)
    cur_k = sum(An[1:1, 1:end])

    writedlm(joinpath("output", "$(t)_$(λ)_$(t)", "$(cur_v)_$(cur_k)_$(t)_$(λ)_$(t).txt"), An, " ")
end
