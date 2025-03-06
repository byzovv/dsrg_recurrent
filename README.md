# Explicit construction of infinite families of strongly regular digraphs

The repository contains a Julia program to generate strongly regular digraphs for the algorithm described in the article "V. A. Byzov, I. A. Pushkarev, Explicit construction of infinite families of strongly regular digraphs with parameters $((v+(2^{n+1}-4)t)2^{n-1}, k+(2^n-2)t, t, \lambda, t)$".

The repository contains:

- directory **input** with the adjacency matrix of the first digraphs in the sequences;
- program **search_dsrg.jl** to search for strongly regular digraphs;
- program **check_dsrg.jl** to check the strongly regular digraphs found;
- archive **output.zip**, containing adjacency matrix of generated digraphs.

---

#### File structure in the catalog input

The catalog **input** contains text files, the names of which have the form of "v_k_t_lambda_mu.txt", where (v, k, t, lambda, mu) are the parameters of the corresponding strongly regular digraphs.

---

#### Program search_dsrg.jl

Before starting the search program for strongly regular digraphs, you need to change in the code the parameters of the first digraph in the sequence:

```julia
first_v = 8
first_k = 4
t = 3
λ = 1
```

The catalog **input** should contain the appropriate text file.

As a result of work, the program creates text files in the catalog **output** in form "v_k_t_lambda_mu.txt".

---

#### Program check_dsrg.jl

The program takes text files in the catalog **output** and checks that the corresponding digraphs are strongly regular.

---

#### Archive output.zip

The archive contains the adjacency matrix of generated digraphs. It is important to note that in the unpacked form, these matrices take up about 8.6 GB of disk space.
