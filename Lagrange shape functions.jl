# This is a generalized function to compute the lagrange shape functions N

using ModelingToolkit

n = 2

# function N(n)
    # N = AbstractArray{n};

    @parameters s

    #=
    for i ∈ 1:n
        for j ∈ 1:n
            if i != j
                N[i] = ((s-j) / (i-j));
            end
        end
    end
    =#

    for i ∈ 1:n
        for j ∈ 1:n
            if i != j
                Ni(s) = (s-j) / (i-j);
            end
        end
    end